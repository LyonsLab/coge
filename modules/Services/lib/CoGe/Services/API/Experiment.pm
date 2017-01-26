package CoGe::Services::API::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(decode_json);
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::Utils;
use CoGe::Core::Experiment qw( delete_experiment );
use CoGe::Core::Favorites;
use CoGe::Core::Metadata qw( create_annotation get_annotation get_annotations );
use CoGe::Services::Auth;
use CoGe::Services::API::Job;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(status => 400, json => { error => { Error => 'Search term is shorter than 3 characters' } });
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

    # Get user's favorites
    my $favorites = CoGe::Core::Favorites->new(user => $user);

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
        restricted => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
        favorited  => $favorites->is_favorite($_) ? Mojo::JSON->true : Mojo::JSON->false
      }
    } @filtered;

    $self->render(json => { experiments => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Invalid input"}
        });
        return;
    }

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $experiment = $self->_get_experiment($id, 0, $db, $user);
    return unless $experiment;

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
            name => $_->name
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

sub fetch_annotations {
    my $self = shift;
    my $id = int($self->stash('id'));
    my ($db) = CoGe::Services::Auth::init($self);

    $self->render(json => get_annotations($id, 'Experiment', $db, 1));
}

sub fetch_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));
    my $aid = int($self->stash('aid'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $experiment = $self->_get_experiment($id, 0, $db, $user);
    return unless $experiment;

    my $annotation = get_annotation($aid, 'Experiment', $db);
    $self->render(json => $annotation) if $annotation;
}

sub add {
    my $self = shift;
    my $data = $self->req->body; #$self->req->json; # mdb replaced 11/22/16 -- req->json hides JSON errors, doing conversion manually prints them to STDERR
    unless ($data) {
        $self->render(status => 400, json => {
            error => { Error => "No request body specified" }
        });
        return;
    }
    $data = decode_json($data);

    # Valid data items
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(status => 400, json => {
            error => { Error => "No data items specified" }
        });
        warn Dumper $data->{source_data};
        return;
    }
    
    # Marshall incoming payload into format expected by Job Submit.
    # Note: This is kind of a kludge -- is there a better way to do this using
    # Mojolicious routing?
    my $request = {
        type => 'load_experiment',
        parameters => $data
    };
    
    return CoGe::Services::API::Job::add($self, $request);
}

sub add_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));

    my $group_name = $self->param('group_name');
    my $image = $self->param('edit_annotation_image');
    my $link = $self->param('link');
    my $text = $self->param('annotation');
    my $type_name = $self->param('type_name');

    my ($db) = CoGe::Services::Auth::init($self);

    create_annotation(
        db => $db,
        group_name => $group_name,
        link => $link,
        target_id => $id,
        target_type => 'experiment',
        text => $text,
        type_name => $type_name
    );
    $self->render(json => { success => Mojo::JSON->true });
}

sub remove {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);

    my $error = delete_experiment($id, $db, $user);
    if ($error) {
        $self->render(status => 400, json => { error => { Error => $error} });
        return;
    }
    $self->render(json => { success => Mojo::JSON->true });
}

sub update {
	my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Invalid input"}
        });
        return;
    }

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $experiment = $self->_get_experiment($id, 1, $db, $user);
    return unless $experiment;

    my $data = $self->req->json;
    if (exists($data->{metadata}->{id})) {
	    delete $data->{metadata}->{id};
    }
	$experiment->update($data->{metadata});
	$self->render(json => { success => Mojo::JSON->true });
}

sub _get_experiment {
    my ($self, $id, $own_or_edit, $db, $user) = @_;
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(status => 404, json => { error => { Error => "Resource not found" } });
        return;
    }
    if ($own_or_edit) {
        unless ($user) {
            $self->render(status => 404, json => { error => { Error => "User not logged in"} });
            return;
        }
        unless ($user->is_owner_editor(experiment => $id)) {
            $self->render(json => { error => { Auth => "Access denied" } }, status => 401);
            return;
        }
    }
    unless ( !$experiment->restricted || (defined $user && $user->has_access_to_experiment($experiment)) ) {
        $self->render(json => { error => { Auth => "Access denied" } }, status => 401);
        return;
    }
    return $experiment;
}

1;
