package CoGe::Services::Data::Group;

use Mojo::Base 'Mojolicious::Controller';
#use IO::Compress::Gzip 'gzip';
use CoGeX;
use CoGe::Accessory::Web;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');
    my $key = $self->param("apiKey");

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);

    # Deny a public group access to groups
    if ($user->is_public) {
        $self->render(json => {
            error => { Auth => "Access denied"}
        }, status => 401);
        return;
    }

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(json => { status => Mojo::JSON->false, error => 'Too many search results'});
        return;
    }

    # Search groups
    my $search_term2 = '%' . $search_term . '%';
    my @groups = $db->resultset("UserGroup")->search(
        \[
            'user_group_id= ? OR name LIKE ? OR description LIKE ?',
            [ 'user_group_id', $search_term ],
            [ 'name', $search_term2 ],
            [ 'description', $search_term2 ],
        ]
    );

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
        role => {
            name => $_->role->name,
            description => $_->role->description
        }
      }
    } @groups;

    $self->render(json => { groups => \@result });
}

sub fetch {
    my $self = shift;
    $self->render(json => { success => Mojo::JSON->true });
}

1;
