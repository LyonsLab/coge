package CoGe::Services::Auth;

use Mojo::UserAgent;
use Data::Dumper;
use URI::Escape::JavaScript qw(unescape);
use JSON qw(decode_json);
use CoGe::Accessory::Web qw(get_defaults add_user parse_proxy_response);

sub init {
    my $self = shift;
    return unless $self;
    #print STDERR Dumper $self, "\n";
    my $username  = $self->param('username');
    my $token     = $self->param('token');
    my $use_cas   = $self->param('use_cas'); # mdb added 7/20/15 for DE
    my $remote_ip = $ENV{REMOTE_ADDR}; #$self->req->env->{HTTP_X_FORWARDED_FOR};
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
    my $user;
    if ($username) {
        $user = $db->resultset('User')->find( { user_name => $username } );
    }

    # Check for existing user session (browser-only)
    my $session_id = unescape($self->cookie($conf->{COOKIE_NAME}));
    if ($session_id) {
        #print STDERR "session_id: ", $session_id, "\n";
        $session_id =~ s/session&//;
        my $session = $db->resultset('UserSession')->find( { session => $session_id } );
        if ($session && $user && $session->user_id == $user->id) {
            print STDERR "CoGe::Services::Auth::init using existing session\n";
            return ( $db, $user, $conf );
        }
    }

    # Otherwise, try to validate user token
    if ($token) {
        my $this_url = $self->req->url->to_abs; # get URL of this request
        my ($uname, $fname, $lname, $email) = validate($username, $token, $use_cas, $this_url);
        unless ($uname) {
            print STDERR "CoGe::Services::Auth::init: token validation failed\n";
            return ( $db, undef, $conf );
        }
    
        # Add new user to DB
        if (!$user) {
            # Add user to database
            print STDERR "CoGe::Services::Auth::init: adding user '", $username, "'\n";
            $user = add_user($coge, $uname, $fname, $lname, $email);
        }
        
        return ( $db, $user, $conf );
    }
    
    # Return unauthenticated response if no token or existing session
    print STDERR "CoGe::Services::Auth::init finished with no authentication\n";
    return ( $db, undef, $conf );
}

sub validate {
    my ($username, $token, $use_cas, $this_url) = @_;
    return unless ($username and $token);
    print STDERR "CoGe::Services::Auth::validate: username=$username token=$token\n";
    
    my ($uname, $fname, $lname, $email);

    # Note: Mojolicious requires IO::Socket::SSL 1.75, do "cpan upgrade IO::Socket::SSL"
    my $ua = Mojo::UserAgent->new;
    
    # CAS Proxy - mdb added 7/20/15 for DE -------------------------------------
    if ($use_cas) {
        # Get URL for CAS
        my $CAS_URL = get_defaults()->{CAS_URL};
        unless ($CAS_URL) {
            print STDERR "CoGe::Services::Auth::validate: missing CAS_URL\n";
            return;
        }
        
        # Validate proxy ticket and get user credentials
        $this_url =~ s/\?.+$//; # remove query params
        my $url = $CAS_URL.'/proxyValidate?service='.$this_url.'&ticket='.$token;
        my $res = $ua->get($url)->res;
        print STDERR Dumper $res, "\n";
        
        ($uname, $fname, $lname, $email) = parse_proxy_response($res->{content}{asset}{content});
        unless ($uname) {
            print STDERR 'CoGe::Services::Auth::validate: CAS failed to authenticate, message=',
                ' url=', $url, "\n";
            return;
        }
    }
    # Agave API ---------------------------------------------------------------
    else {
        # Get URL for Agave User API endpoint
        my $USER_API_URL = get_defaults()->{USER_API_URL};
        unless ($USER_API_URL) {
            print STDERR "CoGe::Services::Auth::validate: missing USER_API_URL\n";
            return;
        }
        
        # Validate token and get user credentials.  We lookup the 'me' profile 
        # for the given token to verify that it belongs to given username.
        # See "Finding yourself" at http://preview.agaveapi.co/documentation/beginners-guides/user-discovery/
        my $url = $USER_API_URL . '/me';
        my $res = $ua->get($url, { Authorization => "Bearer $token" })->res;
        #print STDERR Dumper $res, "\n";
        unless ($res and $res->{message} eq 'OK') {
            print STDERR 'CoGe::Services::Auth::validate: user agent error, message=',
                ($res ? $res->{message} : 'undef'),
                ' url=', $url, "\n";
            print STDERR Dumper $res, "\n" if ($res);
            return;
        }
        
        # Extract user information and verify that the given username owns the given token
        my $authResponse = decode_json('{"' . $res->{content}->{post_buffer}); #FIXME this is a hack because the response is malformed for some unknown reason
        #print STDERR Dumper $authResponse, "\n";
        unless ($authResponse && $authResponse->{status} =~ /success/i &&
                $authResponse->{result} && $authResponse->{result}->{username} eq $username)
        {
            print STDERR 'CoGe::Services::Auth::validate: Agave failed to authenticate, message=',
                ($authResponse ? $authResponse->{message} : 'undef'),
                ' username=',
                ($authResponse ? $authResponse->{result}->{username} : 'undef'),
                ' url=', $url, "\n";
            print STDERR Dumper $authResponse, "\n" if ($authResponse);
            return;
        }
        
        $uname = $authResponse->{result}->{username};
        $fname = $authResponse->{result}->{firstName};
        $lname = $authResponse->{result}->{lastName};
        $email = $authResponse->{result}->{email};
    }

    print STDERR "CoGe::Services::Auth::validate: success! " . ($uname ? $uname : '') . ' ' . ($this_url ? $this_url : '') . "\n";
    return ($uname, $fname, $lname, $email);
}

1;
