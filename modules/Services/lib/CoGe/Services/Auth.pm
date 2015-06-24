package CoGe::Services::Auth;

use Mojo::UserAgent;
use Data::Dumper;
use CoGe::Accessory::Web;
use URI::Escape::JavaScript qw(unescape);

#our $AUTH_API_URL = 'foundation.iplantc.org/auth-v1'; # mdb removed COGE-581
#our $USER_API_URL = 'user.iplantcollaborative.org/api/v1'; # mdb removed COGE-581 -- moved into coge.conf

sub init {
    my $mojo = shift;
    return unless $mojo;
    my $username = $mojo->param('username');
    my $token    = $mojo->param('token');
    my $remote_ip = $ENV{REMOTE_ADDR}; #$mojo->req->env->{HTTP_X_FORWARDED_FOR};
    print STDERR "CoGe::Services::Auth::init: username=", ($username ? $username : ''), 
                 " token=", ($token ? $token : ''), 
                 " remote_ip=", ($remote_ip ? $remote_ip : ''), "\n";

    # Get config
    my $conf = get_defaults();
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

    # Get user from DB
    my $user = $db->resultset('User')->find( { user_name => $username } );

    # Check for existing user session
    my $session_id = unescape($mojo->cookie($conf->{COOKIE_NAME}));
    $session_id =~ s/session&//;

    my $session = $db->resultset('UserSession')->find( { session => $session_id } );

    # TODO add expiration to session table and check it here

    # Use existing session if possible
    if ($session && $user && $session->user_id == $user->id) {
        print STDERR "CoGe::Services::Auth::init using existing session\n";
        return ( $db, $user, $conf );
    }

    # Otherwise, try to validate user token
    unless (validate($username, $token)) {
        print STDERR "CoGe::Services::Auth::init: validation failed\n";
        return ( $db, undef, $conf );
    }

    # Add user to DB
    if (!$user) {
        # TODO add user to database -- right now the user must login via web site before they can use the API
        print STDERR "CoGe::Services::Auth::init: failed to find user '", $username, "'\n";
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

    # Note: Mojolicious requires IO::Socket::SSL 1.75, do "cpan upgrade IO::Socket::SSL"
    my $ua = Mojo::UserAgent->new;
    
    # Agave API ---------------------------------------------------------------
    my $USER_API_URL = get_defaults()->{USER_API_URL}; # test token against API
    unless ($USER_API_URL) {
        print STDERR "CoGe::Services::Auth::validate: missing USER_API_URL\n";
        return;
    }
    
    my $url = $USER_API_URL . '/' . $username;
    my $res = $ua->get($url, { Authorization => "Bearer $token" })->res;
    #print STDERR Dumper $res, "\n";
    unless ($res and $res->{message} eq 'OK') {
        print STDERR 'CoGe::Services::Auth::validate: user agent error, message=',
            ($res ? $res->{message} : 'undef'),
            ' url=', $url, "\n";
        return;
    }
    my $authResponse = $res->{content}->{post_buffer};
    unless ($authResponse and $authResponse =~ /success/) { #FIXME this is a hack because the response is weird
        print STDERR 'CoGe::Services::Auth::validate: failed to authenticate, message=',
            ($authResponse ? $authResponse->{message} : 'undef'),
            ' url=', $url, "\n";
        return;
    }
    
    # Foundation API # mdb removed COGE-581 -----------------------------------
#    my $url = 'https://'.$username.':'.$token.'@'.$AUTH_API_URL.'/list';
#    my $res = $ua->get($url)->res;
#    #print STDERR Dumper $res, "\n";
#    unless ($res and $res->{message} eq 'OK') {
#        print STDERR 'CoGe::Services::Auth::validate: user agent error, message=',
#            ($res ? $res->{message} : 'undef'),
#            ' url=', $url, "\n";
#        return;
#    }
#    my $authResponse = $res->json;
#    unless ($authResponse and $authResponse->{status} eq 'success') {
#        print STDERR 'CoGe::Services::Auth::validate: failed to authenticate, message=',
#            ($authResponse ? $authResponse->{message} : 'undef'),
#            ' url=', $url, "\n";
#        return;
#    }

    print STDERR "CoGe::Services::Auth::validate: success!\n";
    return 1;
}

1;
