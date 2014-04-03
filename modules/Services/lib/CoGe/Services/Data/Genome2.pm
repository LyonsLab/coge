package CoGe::Services::Data::Genome2;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Accessory::Web;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');
    my $key = $self->param("apiKey");

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(json => { status => Mojo::JSON->false, error => 'Too many search results'});
        return;
    }

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);

    # Search genomes
    my $search_term2 = '%' . $search_term . '%';
    my @genomes= $db->resultset("Genome")->search(
        \[
            'genome_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'genome_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter response
    my @filtered = grep {
        $user->has_access_to_genome($_);
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
    my $key = $self->param("apiKey");

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);

    my $genome = $db->resultset("Genome")->find($id);

    unless (defined $genome) {
        $self->render(json => {
            error => { Error => "Item not found"}
        });
        return;
    }

    unless ($user->has_access_to_genome($genome)) {
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
