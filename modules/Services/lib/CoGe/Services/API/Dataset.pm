package CoGe::Services::API::Dataset;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Services::Error;

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    my $dataset = $db->resultset("Dataset")->find($id);
    unless (defined $dataset) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    unless ( !$dataset->restricted || (defined $user && $user->has_access_to_dataset($dataset)) ) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    $self->render(json => {
        id => int($dataset->id),
        name => $dataset->name,
        description => $dataset->description,
        date => $dataset->date,
        link => $dataset->link,
        version => $dataset->version,
        organism_id  => int($dataset->organism->id),
        data_source => {
            name => $dataset->source->name,
            description => $dataset->source->description,
            link => $dataset->source->link,
        },
        chromosomes => $dataset->chromosome_count(ftid => 4),
        restricted => $dataset->restricted ? Mojo::JSON->true : Mojo::JSON->false,
    });
}

sub genomes {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    #my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);
    my ($db, $user) = CoGe::Services::Auth::init($self);

    my $dataset = $db->resultset("Dataset")->find($id);
    unless (defined $dataset) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    unless ( !$dataset->restricted || (defined $user && $user->has_access_to_dataset($dataset)) ) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Filter response
    my @filtered = grep {
        !$_->restricted || (defined $user && $user->has_access_to_genome($_))
    } $dataset->genomes;

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
      }
    } @filtered;

    return $self->render(json => { genomes => \@result });
}

1;
