package CoGe::Services::Data::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Data::Dumper;
#use IO::Compress::Gzip 'gzip';
use CoGeX;
use CoGe::Accessory::Utils;
use CoGe::Services::Auth;
use CoGe::Services::Data::Job;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Search experiments
    my $search_term2 = '%' . $search_term . '%';
    my @experiments = $db->resultset("Experiment")->search(
        \[
            'experiment_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'experiment_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter response
    my @filtered = grep {
         !$_->restricted || (defined $user && $user->has_access_to_experiment($_))
    } @experiments;

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
      }
    } @filtered;

    $self->render(json => { experiments => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get experiment
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(json => {
            error => { Error => "Item not found" }
        });
        return;
    }

    # Check permissions
    unless ( !$experiment->restricted || (defined $user && $user->has_access_to_experiment($experiment)) ) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group
        }
    } $experiment->annotations;

    # Format types
    my @types = map {
        {
            name => $_->name,
            description => $_->description
        }
    } $experiment->types;

    $self->render(json => {
        id => int($experiment->id),
        name => $experiment->name,
        description => $experiment->description,
        version => $experiment->version,
        genome_id  => int($experiment->genome->id),
        source => {
            name => $experiment->source->name,
            description => $experiment->source->description,
            link => $experiment->source->link
        },
        types => \@types,
        additional_metadata => \@metadata,
        restricted => $experiment->restricted ? Mojo::JSON->true : Mojo::JSON->false,
    });
}

sub add {
    my $self = shift;
    my $data = $self->req->json;
    #print STDERR "CoGe::Services::Data::Experiment::add\n", Dumper $data, "\n";

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required to add experiment
    unless (defined $user) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        });
        return;
    }

    # Valid data items
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(json => {
            error => { Error => "No data items specified" }
        });
        return;
    }
    
    # Marshall incoming payload into format expected by Job Submit.
    # Note: This is kind of a kludge -- is there a better way to do this using
    # Mojolicious routing?
    my $request = {
        type => 'load_experiment',
        parameters => $data
    };
    
    return CoGe::Services::Data::Job::add($self, $request);
}

1;
