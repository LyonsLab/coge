package CoGe::Services::Auth;

use Mojo::UserAgent;
use Data::Dumper;
use URI::Escape::JavaScript qw(unescape);
use JSON qw(decode_json);
use CoGe::Accessory::Web qw(get_defaults add_user parse_proxy_response jwt_decode_token);
use File::Spec::Functions qw(catfile);

sub init {
    my $self = shift;
    return unless $self;
    #print STDERR Dumper $self, "\n";
    my $username  = $self->param('username');
    my $token     = $self->param('token');
    my $token2    = $self->req->headers->header('x-iplant-de-jwt'); # mdb added 9/23/15 for DE
    my $remote_ip = $ENV{REMOTE_ADDR}; #$self->req->env->{HTTP_X_FORWARDED_FOR};
#    warn "CoGe::Services::Auth::init";
#    warn "username=" . ($username ? $username : '');
#    warn "token=" . ($token ? $token : '');
#    warn "token2=" . ($token2 ? $token2 : '');
#    warn "remote_ip=" . ($remote_ip ? $remote_ip : '');

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

    # Check for existing user session (cookie enabled browser only)
    print STDERR Dumper $self->cookie($conf->{COOKIE_NAME});
    my $session_id = unescape($self->cookie($conf->{COOKIE_NAME}));
    if ($session_id) {
        print STDERR "session_id: ", $session_id, "\n";
        $session_id =~ s/session&//;
        my $session = $db->resultset('UserSession')->find( { session => $session_id } );
        if ($session) {# && $user && $session->user_id == $user->id) { # mdb changed 3/7/16 for hypnotoad
            $user = $db->resultset('User')->find($session->user_id); # mdb added 3/7/16 for hypnotoad
            print STDERR "CoGe::Services::Auth::init using existing session for user '", $user->name, "'\n";
            return ( $db, $user, $conf );
        }
    }

    # Otherwise, try to validate user token
    if ($token || $token2) {
        my ($uname, $fname, $lname, $email);
        if ($token) { # Agave
            ($uname, $fname, $lname, $email) = validate_agave($username, $token);
        }
        else { # DE JWT
            ($uname, $fname, $lname, $email) = validate_jwt($token2, $conf);
        }
        
        unless ($uname) {
            print STDERR "CoGe::Services::Auth::init: token validation failed\n";
            return ( $db, undef, $conf );
        }
    
        # Add new user to DB
        if (!$user) {
            print STDERR "CoGe::Services::Auth::init: adding user '", $uname, "'\n";
            $user = add_user($db, $uname, $fname, $lname, $email);
        }
        
        return ( $db, $user, $conf );
    }
    
    # Return unauthenticated response if no token or existing session
    print STDERR "CoGe::Services::Auth::init finished with no authentication\n";
    return ( $db, undef, $conf );
}

sub validate_jwt {
    my $token = shift;
    my $conf = shift;
    return unless $token;
    print STDERR "CoGe::Services::Auth::validate_jwt\n";
    
    # Get path to DE public key file
    my $de_public_key_path = catfile($conf->{RESOURCEDIR}, $conf->{DE_PUBLIC_KEY});
    unless ($de_public_key_path) {
        print STDERR "CoGe::Services::Auth::init: missing DE_PUBLIC_KEY in config file\n";
        return;
    }
    
    # Decode token and get payload
    my $claims = jwt_decode_token($token, $de_public_key_path);
    unless ($claims) {
        print STDERR "CoGe::Services::Auth::validate_jwt: JWT token decoding failed\n";
        return;
    }
        
    my $uname = $claims->{'sub'};
    my $fname = $claims->{'given_name'};
    my $lname = $claims->{'family_name'};
    my $email = $claims->{'email'};
        
    print STDERR "CoGe::Services::Auth::validate_jwt: success! ", ($uname ? $uname : ''), "\n";
    return ($uname, $fname, $lname, $email);
}

sub validate_agave {
    my ($username, $token) = @_;
    return unless ($username and $token);
    print STDERR "CoGe::Services::Auth::validate_agave: username=$username token=$token\n";
    
    my ($uname, $fname, $lname, $email);

    # Note: Mojolicious requires IO::Socket::SSL 1.75, do "cpan upgrade IO::Socket::SSL"
    my $ua = Mojo::UserAgent->new;

    # CAS Proxy - mdb added 7/20/15 for DE -------------------------------------
#    if ($token_type eq 'cas') {
#        # Get URL for CAS
#        my $CAS_URL = get_defaults()->{CAS_URL};
#        unless ($CAS_URL) {
#            print STDERR "CoGe::Services::Auth::validate: missing CAS_URL\n";
#            return;
#        }
#        
#        # Validate proxy ticket and get user credentials
#        $this_url =~ s/\?.+$//; # remove query params
#        my $url = $CAS_URL.'/proxyValidate?service='.$this_url.'&ticket='.$token;
#        my $res = $ua->get($url)->res;
#        print STDERR Dumper $res, "\n";
#        
#        ($uname, $fname, $lname, $email) = parse_proxy_response($res->{content}{asset}{content});
#        unless ($uname) {
#            print STDERR 'CoGe::Services::Auth::validate_agave: CAS failed to authenticate, message=',
#                ' url=', $url, "\n";
#            return;
#        }
#    }

    # Agave API (default) ------------------------------------------------------
    
    # Get URL for Agave User API endpoint
    my $USER_API_URL = get_defaults()->{USER_API_URL};
    unless ($USER_API_URL) {
        print STDERR "CoGe::Services::Auth::validate_agave: missing USER_API_URL\n";
        return;
    }

    # Validate token and get user credentials.  We lookup the 'me' profile 
    # for the given token to verify that it belongs to given username.
    # See "Finding yourself" at http://preview.agaveapi.co/documentation/beginners-guides/user-discovery/
    my $url = $USER_API_URL . '/me';
    my $res = $ua->get($url, { Authorization => "Bearer $token" })->res;
    unless ($res and $res->{message} eq 'OK') {
        print STDERR 'CoGe::Services::Auth::validate: user agent error, message=',
            ($res ? $res->{message} : 'undef'),
            ' url=', $url, "\n";
        print STDERR Dumper $res, "\n" if ($res);
        return;
    }
    
    # Extract user information and verify that the given username owns the given token
    #print STDERR Dumper $res->body, "\n";
    my $authResponse = decode_json($res->body);
    unless ($authResponse && $authResponse->{status} =~ /success/i &&
            $authResponse->{result} && $authResponse->{result}->{username} eq $username)
    {
        print STDERR 'CoGe::Services::Auth::validate_agave: Agave failed to authenticate, message=',
            ($authResponse ? $authResponse->{message} : 'undef'),
            ' url=', $url, "\n";
        print STDERR Dumper $authResponse, "\n" if ($authResponse);
        return;
    }
    
    $uname = $authResponse->{result}->{username};
    $fname = $authResponse->{result}->{firstName};
    $lname = $authResponse->{result}->{lastName};
    $email = $authResponse->{result}->{email};

    print STDERR "CoGe::Services::Auth::validate_agave: success! ", ($uname ? $uname : ''), "\n";
    return ($uname, $fname, $lname, $email);
}

1;
