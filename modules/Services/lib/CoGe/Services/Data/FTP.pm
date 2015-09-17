package CoGe::Services::Data::FTP;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Accessory::Web qw(ftp_get_path);
use CoGe::Services::Auth;
use Data::Dumper;

sub list {
    my $self = shift;
    #my $url = $self->stash('url');
    my $url = $self->param('url');
    print STDERR "FTP::list ", $url, "\n";

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(json => { error => { Error => 'Access denied' } });
        return;
    }

    # Fetch directory listing
    my $files = ftp_get_path(url => $url);
    unless ($files) {
        $self->render(json => { error => { FTP => 'Not found' } });
        return;
    }
    
    $self->render(json => { url => $url, items => $files });
}

1;
