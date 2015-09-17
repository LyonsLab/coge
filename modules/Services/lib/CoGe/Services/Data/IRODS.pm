package CoGe::Services::Data::IRODS;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Accessory::IRODS;
use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_irods_path get_irods_file);
use CoGe::Services::Auth;
use Data::Dumper;

sub list {
    my $self = shift;
    my $path = $self->stash('path');
    #print STDERR "IRODS::list ", $path, "\n";

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(json => { error => { Error => 'Access denied' } });
        return;
    }

    # Fetch directory listing
    my $result = get_irods_path($path, $user->name);
    unless ($result) {
        $self->render(json => { error => { IRODS => 'Access denied' } });
        return;
    }
    
    my $error  = $result->{error};
    if ($error) {
        $self->render(json => { error => { IRODS => $error } });
        return;
    }

    $self->render(json => { path => $result->{path}, items => $result->{items} });
}

# mdb removed 8/24/15 -- not used
#sub fetch {
#    my $self = shift;
#    my $path = $self->stash('path');
#    my $load_id = $self->param('load_id');
#    $load_id = get_unique_id() unless $load_id;
#    #print STDERR "IRODS::fetch ", $path, "\n";
#    
#    # Authenticate user and connect to the database
#    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
#    unless ($user) {
#        $self->render(json => { error => { Error => 'Access denied' } });
#        return;
#    }
#    
#    $path = unescape($path);
#    my $uploadpath = get_upload_path($user->name, $load_id);
#    my $result = get_irods_file($path, $uploadpath);
#
#    $self->render(json => $result );
#}

1;
