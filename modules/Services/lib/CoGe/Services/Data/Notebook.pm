package CoGe::Services::Data::Notebook;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
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
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Search genomes
    my $search_term2 = '%' . $search_term . '%';
    my @notebooks = $db->resultset("List")->search(
        \[
            'list_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'list_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter response
    my @filtered = grep {
        !$_->restricted || (defined $user && $user->has_access_to_list($_))
    } @notebooks;

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description
      }
    } @filtered;

    $self->render(json => { notebooks => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(json => {
            error => { Error => "Item not found"}
        });
        return;
    }

    unless ( !$notebook->restricted || (defined $user && $user->has_access_to_genome($notebook)) ) {
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
    } $notebook->annotations;

    # Format items
    my @items = map {
        { type => 'genome', id => $_->id }
    } $notebook->genomes;

    push @items, map {
        { type => 'experiment', id => $_->id }
    } $notebook->experiments;

    # Format response
    $self->render(json => {
        id => int($notebook->id),
        name => $notebook->name,
        description => $notebook->description,
        type => $notebook->type->name,
        restricted => $notebook->restricted ? Mojo::JSON->true : Mojo::JSON->false,
        metadata => \@metadata,
        items => \@items
    });
}

1;
