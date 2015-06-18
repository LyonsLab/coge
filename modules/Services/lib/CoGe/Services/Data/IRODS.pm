package CoGe::Services::Data::IRODS;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Accessory::IRODS;
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

    # Setup path
    #TODO set path to home if not specified
    #my $username = $user->name;
    #my $basepath = $conf->{IRODSDIR};
    #$basepath =~ s/\<USER\>/$username/;
    #$path = $basepath unless $path;
    $path = '/' . $path;

    # Fetch directory listing
    my $result = CoGe::Accessory::IRODS::irods_ils($path);
    my $error  = $result->{error};
    if ($error) {
        $self->render(json => { error => { IRODS => $error } });
        return;
    }

    $self->render(json => { path => $path, items => $result->{items} });
}

sub fetch {
    my $self = shift;
    my $path = $self->stash('path');
    #print STDERR "IRODS::fetch ", $path, "\n";
    
    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(json => { error => { Error => 'Access denied' } });
        return;
    }
    
    $path = unescape($path);
    my ($filename)   = $path =~ /([^\/]+)\s*$/;
    my ($remotepath) = $path =~ /(.*)$filename$/;

    my $localpath     = catdir('irods', $remotepath);
    my $localfullpath = catdir($TEMPDIR, $localpath);
    $localpath = catfile($localpath, $filename);
    my $localfilepath = catfile($localfullpath, $filename);

    my $do_get = 1;

    #   if (-e $localfilepath) {
    #       my $remote_chksum = irods_chksum($path);
    #       my $local_chksum = md5sum($localfilepath);
    #       $do_get = 0 if ($remote_chksum eq $local_chksum);
    #       print STDERR "$remote_chksum $local_chksum\n";
    #   }

    if ($do_get) {
        mkpath($localfullpath);
        CoGe::Accessory::IRODS::irods_iget( $path, $localfullpath );
    }

    $self->render(json => { path => $localpath, size => -s $localfilepath } );
}

1;
