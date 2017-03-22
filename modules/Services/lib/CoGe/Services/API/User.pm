package CoGe::Services::API::User;

use Mojo::Base 'Mojolicious::Controller';
#use IO::Compress::Gzip 'gzip';
use Data::Dumper;
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Services::Error;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Authenticate user and connect to the database
    #my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Deny a public user access to users
    unless ($user && !$user->is_public) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Search users
    my $search_term2 = '%' . $search_term . '%';
    my @users = $db->resultset("User")->search(
        \[
            'user_id = ? OR user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
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
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    #my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $key);
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Deny a public user access to any user
    unless ($user && !$user->is_public) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Get requested user
    my $fetched_user = $db->resultset("User")->find($id);
    unless (defined $fetched_user) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Restrict a user's access to themselves
    if ($user->id ne $fetched_user->id) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    $self->render(json => {
        id => int($fetched_user->id),
        user_name => $fetched_user->user_name,
        first_name => $fetched_user->first_name,
        last_name => $fetched_user->last_name,
        #email => $fetched_user->email,
        description => $fetched_user->description,
    });
}

sub items {
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

    # Deny a public user access to any user
    unless ($user && !$user->is_public) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    my $fetched_user = $db->resultset("User")->find($id);
    unless (defined $fetched_user) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Restrict a user's access to themselves
    if ($user->id ne $fetched_user->id) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Format experiments
    my @experiments = map {
        {
            id => int($_->id),
            type => "experiment",
        }
    } $fetched_user->experiments;

    # Format Genomes
    my @genomes = map {
        {
            id => int($_->id),
            type => "genome",
        }
    } $fetched_user->genomes;

    # Format Lists
    my @notebooks = map {
        {
            id => int($_->id),
            type => "notebook"
        }
    } $fetched_user->lists;

    $self->render(json => {
        items => [@experiments, @genomes, @notebooks]
    });
}

1;
