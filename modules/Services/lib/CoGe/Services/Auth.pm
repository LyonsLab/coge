package CoGe::Services::Auth;

use Mojo::UserAgent;
use Data::Dumper;
use CoGe::Accessory::Web;

# FIXME move these into coge.conf someday
our $AUTH_API_URL = 'foundation.iplantc.org/auth-v1';
our $USER_API_URL = 'user.iplantcollaborative.org/api/v1';

sub init {
    my $mojo = shift;
    return unless $mojo;
    my $opts = shift;
    my $username = $mojo->param('username');
    my $token    = $mojo->param('token');
    my $remote_ip = $mojo->req->env->{HTTP_X_FORWARDED_FOR};
    print STDERR "CoGe::Services::Auth::init: username=", ($username ? $username : ''), " token=", ($token ? $token : ''), " remote_ip=", ($remote_ip ? $remote_ip : ''), "\n";

    # Get config
    my $conf = CoGe::Accessory::Web::get_defaults();
    unless (defined $conf) {
        print STDERR "CoGe::Services::Auth::init: couldn't load config file\n";
        return;
    }

    # Connect to DB
    my $db = CoGeX->dbconnect($conf);
#    if ($debug) { # enable ORM debugging if requested
#        $db->storage->debugobj(new DBIxProfiler());
#        $db->storage->debug(1);
#    }
    unless (defined $db) {
        print STDERR "CoGe::Services::Auth::init: couldn't connect to database\n";
        return;
    }

    # Check for existing user session
    #my $session_id = CoGe::Accessory::Web::get_session_id($username, $remote_ip);
    my $session;# = $db->resultset('UserSession')->find( { session => $session_id } );
    # TODO add expiration to session table and check it here

    # Otherwise, validate user
    unless ($session or validate($username, $token)) {
        print STDERR "CoGe::Services::Auth::init: validation failed\n";
        return ( $db, undef, $conf );
    }

    # Get user from DB
    my $user = $db->resultset('User')->find( { user_name => $username } );

    # Otherwise, get user info and add to DB
    if (!$user) {
        # TODO get user info from iPlant user API and add user to database
        return ( $db, undef, $conf ); # tempfix
    }

    # Lastly, create or renew session
    # TODO

    return ( $db, $user, $conf );
}

sub validate {
    my ($username, $token) = @_;
    return unless ($username and $token);
    print STDERR "CoGe::Services::Auth::validate: username=$username token=$token\n";

    # Validate ticket
    # Note: Mojolicious requires IO::Socket::SSL 1.75, do "cpan upgrade IO::Socket::SSL"
    my $ua = Mojo::UserAgent->new;
    my $url = 'https://'.$username.':'.$token.'@'.$AUTH_API_URL.'/list';
    my $res = $ua->get($url)->res;
    #print STDERR Dumper $res, "\n";
    unless ($res and $res->{message} eq 'OK') {
        print STDERR 'CoGe::Services::Auth::validate: user agent error, message=',
            ($res ? $res->{message} : 'undef'),
            ' url=', $url, "\n";
        return;
    }
    my $authResponse = $res->json;
    unless ($authResponse and $authResponse->{status} eq 'success') {
        print STDERR 'CoGe::Services::Auth::validate: failed to authenticate, message=',
            ($authResponse ? $authResponse->{message} : 'undef'),
            ' url=', $url, "\n";
        return;
    }

    return 1;
}

sub get_user {
    my ($username) = @_;

    # Get user information if not already in our database
    $url = 'https://coge:150b316673023dff56f7ba125df21c4f@'.$USER_API_URL.'/users/username/'.$username;
    $res = $ua->get($url)->res;
    print STDERR Dumper $res, "\n";
    unless ($res and $res->{message} eq 'OK') {
        print STDERR 'CoGe::Services::Auth::validate: user agent error, message=',
            ($res ? $res->{message} : 'undef'),
            ' url=', $url, "\n";
        return;
    }
    my $userResponse = $res->json;
    unless ($userResponse) {
        print STDERR 'CoGe::Services::Auth::validate: failed to authenticate, message=',
            ($userResponse ? $userResponse->{message} : 'undef'),
            ' url=', $url, "\n";
        return;
    }
}

1;
