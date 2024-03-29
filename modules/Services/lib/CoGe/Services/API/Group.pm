package CoGe::Services::API::Group;

use Mojo::Base 'Mojolicious::Controller';
use Data::Dumper;
#use IO::Compress::Gzip 'gzip';

use CoGeX;
use CoGe::Services::Auth;
use CoGe::Services::Error;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Deny a public user access to groups
    unless ($user && !$user->is_public) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(API_STATUS_SEARCHTERM);
        return;
    }

    # Search groups
    my $search_term2 = '%' . $search_term . '%';
    my @groups = $db->resultset("UserGroup")->search(
        \[
            'user_group_id = ? OR name LIKE ? OR description LIKE ?',
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
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $group = $self->_get_group($id, 1, $db, $user);
    return unless $group;

    $self->render(json => {
        id => int($group->id),
        name => $group->name,
        description => $group->description,
        role => {
            name => $group->role->name,
            description => $group->role->description
        },
        users => [ map { int($_->id) } $group->users ]
    });
}

sub items {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $group = $self->_get_group($id, 1, $db, $user);
    return unless $group;

    # Restrict access to users in the group
    unless ($group->has_member($user)) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Format experiments
    my @experiments = map {
        {
            id => int($_->id),
            type => "experiment",
        }
    } $group->experiments;

    # Format Genomes
    my @genomes = map {
        {
            id => int($_->id),
            type => "genome",
        }
    } $group->genomes;

    # Format Lists
    my @notebooks = map {
        {
            id => int($_->id),
            type => "notebook"
        }
    } $group->lists;

    $self->render(json => {
        items => [@experiments, @genomes, @notebooks]
    });
}

sub remove_users {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $group = $self->_get_group($id, 1, $db, $user);
    return unless $group;

    if ($group->locked && !$user->is_admin) {
        $self->render(API_STATUS__BAD_REQUEST("Group locked"));
        return;
    }

    my $data = $self->req->json;
    my @users;
    foreach my $user_id (@{$data->{users}}) {
        my $uc = $db->resultset('UserConnector')->find({
            parent_id   => $user_id,
            parent_type => 5,
            child_id    => $id,
            child_type  => 6
        });
        if ($uc && !$uc->role->is_owner) {
            $uc->delete;
            my $u = $db->resultset('User')->find($user_id);
            CoGe::Accessory::Web::log_history(
                db          => $db,
                user_id     => $user->id,
                description => 'removed user ' . $u->info . ' from group ' . $group->info_html,
                page        => 'API',
                parent_id   => $id,
                parent_type => 6 #FIXME magic number
            );
        }
    }
    $self->render(json => { success => Mojo::JSON->true });
}

sub _get_group {
    my ($self, $id, $own_or_edit, $db, $user) = @_;

    # Deny a public user access to any group
    unless ($user && !$user->is_public) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    my $group = $db->resultset("UserGroup")->find($id);
    unless (defined $group) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    if ($own_or_edit) {
        unless ($group->is_editable($user)) {
            $self->render(API_STATUS_UNAUTHORIZED);
            return;
        }
    }
    return $group;
}

1;
