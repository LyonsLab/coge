package CoGe::Services::Data::Notebook;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Core::Notebook;

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

    # Search notebooks and filter based on user permissions
    my $aNotebooks = search_notebooks(db => $db, search_term => $search_term, user => $user);

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description
      }
    } @$aNotebooks;

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
