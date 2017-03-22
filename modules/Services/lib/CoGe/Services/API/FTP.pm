package CoGe::Services::API::FTP;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Accessory::Web qw(ftp_get_path);
use CoGe::Services::Auth;
use CoGe::Services::Error;

sub list {
    my $self = shift;
    my $url = $self->param('url');
    my $dirs = $self->param('dirs');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Fetch directory listing
    my $files = ftp_get_path(url => $url, dirs => $dirs);
    unless ($files) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }
    
    $self->render(json => { url => $url, items => $files });
}

1;
