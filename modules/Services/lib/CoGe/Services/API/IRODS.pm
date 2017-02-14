package CoGe::Services::API::IRODS;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Accessory::IRODS;
use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_irods_path get_irods_file irods_mkdir irods_rm);
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

sub mkdir {
    my $self = shift;
    my $path = $self->param('path');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(json => { error => { Error => 'Access denied' } });
        return;
    }

    my $error = irods_mkdir($path);
    if ($error) {
        $self->render(json => { error => { Error => $error } });
        return;
    }
    $self->render(json => { success => Mojo::JSON->true });
}

1;
