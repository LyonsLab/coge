package CoGe::Services::Data::User;

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

    # Deny a public user access to users
    if ($user->is_public) {
        $self->render(json => {
            error => { Auth => "Access denied"}
        }, status => 401);
        return;
    }

    # Search users
    my $search_term2 = '%' . $search_term . '%';
    my @users = $db->resultset("User")->search(
        \[
            'user_id= ? OR user_name LIKE ? OR first_name LIKE ?
                OR last_name LIKE ?',
            [ 'user_id', $search_term ],
            [ 'user_name', $search_term2 ],
            [ 'first_name', $search_term2 ],
            [ 'last_name', $search_term2 ]
        ]
    );

    # Format response
    my @result = map {
      {
        id => int($_->id),
        user_name => $_->user_name,
        first_name => $_->first_name,
        last_name => $_->last_name,
        description => $_->description
      }
    } @users;

    $self->render(json => { users => \@result });
}

sub fetch {
    my $self = shift;
    $self->render(json => { success => Mojo::JSON->true });
}

1;
