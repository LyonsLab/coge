package CoGe::Services::Data::Genome2;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Core::Genome qw(search_genomes);
use CoGe::Services::Auth;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    #my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);
    my ($db, $user) = CoGe::Services::Auth::init($self);

    my @genomes = search_genomes(db => $db, user => $user, search => $search_term);

    # Filter response
    my @filtered = grep {
        !$_->restricted || (defined $user && $user->has_access_to_genome($_))
    } @genomes;

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
      }
    } @filtered;

    $self->render(json => { genomes => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    #my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);
    my ($db, $user) = CoGe::Services::Auth::init($self);

    my $genome = search_genomes(db => $db, gid => $id);

    unless (defined $genome) {
        $self->render(json => {
            error => { Error => "Item not found"}
        });
        return;
    }

    unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
        $self->render(json => {
            error => { Auth => "Access denied"}
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
    } $genome->annotations;

    $self->render(json => {
        id => int($genome->id),
        name => $genome->name,
        description => $genome->description,
        link => $genome->link,
        version => $genome->version,
        organism_id  => int($genome->organism->id),
        sequence_type => {
            name => $genome->type->name,
            description => $genome->type->description,
        },
        experiments => [ map { int($_->id) } $genome->experiments ],
        metadata => \@metadata,
        restricted => $genome->restricted ? Mojo::JSON->true : Mojo::JSON->false,
    });
}

1;
