package CoGe::Accessory::Web;

use strict;
use CoGeX;
use Data::Dumper;
use base 'Class::Accessor';
use CGI::Carp('fatalsToBrowser');
use CGI;
use CGI::Cookie;
use DBIxProfiler;
use File::Basename;
use File::Temp;
use LWP::Simple qw(!get !head !getprint !getstore !mirror);
use LWP::UserAgent;
use HTTP::Request;
use XML::Simple;
use CoGe::Accessory::LogUser;
use Digest::MD5 qw(md5_base64);
use POSIX;
use IPC::System::Simple qw(capture system $EXITVAL EXIT_ANY);
use Mail::Mailer;

BEGIN {
	use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK $Q $cogex $TEMPDIR $BASEDIR);
	require Exporter;

	$BASEDIR = "/opt/apache/CoGe/";
	$VERSION = 0.1;
	$TEMPDIR = $BASEDIR . "tmp";
	@ISA     = ( @ISA, qw (Exporter) );

	#Give a hoot don't pollute, do not export more than needed by default
	@EXPORT = qw ()
	  ; #qw (login write_log read_log check_taint check_filename_taint save_settings load_settings reset_settings initialize_basefile);

#    $cogex = CoGeX->dbconnect();
#    $cogex->storage->debugobj(new DBIxProfiler());
#    $cogex->storage->debug(1);
	__PACKAGE__->mk_accessors('restricted_orgs', 'basefilename',  'basefile', 'logfile',  'sqlitefile');
}

sub get_defaults {
	my ( $self, $param_file ) = self_or_default(@_);
	$param_file = $BASEDIR . "coge.conf" unless defined $param_file;
	unless ( -r $param_file ) {
		print STDERR
qq{Either no parameter file specified or unable to read paramer file ($param_file).
A valid parameter file must be specified or very little will work!};
		return 0;
	}
	open( IN, $param_file );
	my %items;
	while (<IN>) {
		chomp;
		next if /^#/;
		next unless $_;
		my ( $name, $path ) = split( /\s+/, $_, 2 );
		$items{$name} = $path;
	}
	close IN;
	return \%items;
}

sub dataset_search_for_feat_name {
	my ( $self, $accn, $num, $dsid, $featid, $coge ) = self_or_default(@_);
	$num = 1 unless $num;
	return (qq{<input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">}, $num)
	  unless $accn;
	my $html;
	my %sources;
	my %restricted_orgs = %{ $self->restricted_orgs } if $self->restricted_orgs;
	my $rs = $coge->resultset('Dataset')->search(
		{ 'feature_names.name' => $accn },
		{
			'join'     => { 'features'    => 'feature_names' },
			'prefetch' => [ 'datasource', 'organism' ]
		}
	);
	while ( my $ds = $rs->next() ) {
		my $name    = $ds->name;
		my $ver     = $ds->version;
		my $desc    = $ds->description;
		my $sname   = $ds->datasource->name;
		my $ds_name = $ds->name;
		my $org     = $ds->organism->name;
		my $title   = "$org: $ds_name ($sname, v$ver)";
		next if $restricted_orgs{$org};
		$sources{ $ds->id } = { title => $title, version => $ver };
	}
	if ( keys %sources ) {
		$html .= qq{
 <SELECT name = "dsid$num" id= "dsid$num" onChange="feat_search(['accn$num','dsid$num', 'args__$num'],['feat$num']);" >
 };
		foreach
		  my $id ( sort { $sources{$b}{version} <=> $sources{$a}{version} }
			keys %sources )
		{
			my $val = $sources{$id}{title};
			$html .= qq{  <option value="$id"};
			$html .= qq{ selected } if $dsid && $id == $dsid;
			$html .= qq{>$val\n};
		}
		$html .= qq{</SELECT>\n};
		my $count = scalar keys %sources;
		$html .= qq{<font class=small>($count)</font>};
	}
	else {
		$html .=
qq{Accession not found <input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">\n};
	}
	return ( $html, $num );
}

sub feat_search_for_feat_name {
	my ( $self, $accn, $dsid, $num, $coge ) = self_or_default(@_);
	return qq{<input type="hidden" id="featid$num">\n} unless $dsid;
	my @feats;
	my $rs = $coge->resultset('Feature')->search(
		{
			'feature_names.name' => $accn,
			'dataset.dataset_id' => "$dsid",
		},
		{
			'join'     => [ 'feature_type', 'dataset', 'feature_names' ],
			'prefetch' => [ 'feature_type', 'dataset' ],
		}
	);
	my %seen;
	while ( my $f = $rs->next() ) {
		next unless $f->dataset->id == $dsid;

		#	next if $f->feature_type->name =~ /CDS/i;
		#	next if $f->feature_type->name =~ /RNA/i;
		push @feats, $f unless $seen{ $f->id };
		$seen{ $f->id } = 1;
	}
	my $html;
	if (@feats) {
		$html .= qq{<SELECT name = "featid$num" id = "featid$num" >};
		foreach my $feat ( sort { $a->type->name cmp $b->type->name } @feats ) {
			my $loc = "(" . $feat->type->name . ") Chr:" . $feat->locations->next->chromosome . " " . $feat->start . "-" . $feat->stop;

			#working here, need to implement genbank_location_string before I can progress.  Need
			$loc =~ s/(complement)|(join)//g;
			my $fid = $feat->id;
			$html .= qq {  <option value="$fid">$loc \n};
		}
		$html .= qq{</SELECT>\n};
		my $count = scalar @feats;
		$html .= qq{<font class=small>($count)</font>};
	}
	else {
		$html .= qq{<input type="hidden" id="featid$num">\n};
	}
	return $html;
}

sub self_or_default {    #from CGI.pm
	return @_
	  if defined( $_[0] )
	  && ( !ref( $_[0] ) )
	  && ( $_[0] eq 'CoGe::Accessory::Web' );
	unless (
		defined( $_[0] )
		&& ( ref( $_[0] ) eq 'CoGe::Accessory::Web'
			|| UNIVERSAL::isa( $_[0], 'CoGe::Accessory::Web' )
		)                # slightly optimized for common case
	  )
	{
		$Q = CoGe::Accessory::Web->new unless defined($Q);
		unshift( @_, $Q );
	}
	return wantarray ? @_ : $Q;
}

sub logout_cas {
	my $self        = shift;
	my %opts        = @_;
	my $cookie_name = $opts{cookie_name};
	my $coge        = $opts{coge};
	my $user        = $opts{user};
	my $form        = $opts{form};          #CGI form for calling page
	my $url         = $opts{this_url};
	my $cas_url = 'http://coge.iplantcollaborative.org/coge/';
	my %cookies = fetch CGI::Cookie;
	$url = $form->url() unless $url;
	my $session = md5_base64( $user->user_name . $ENV{REMOTE_ADDR} );
	$session =~ s/\+/1/g;
	($session) = $coge->resultset('UserSession')->find( { session => $session } );
	$session->delete if $session;
	print "Location: " . $form->redirect("https://auth.iplantcollaborative.org/cas/logout?service=" . $url . "&gateway=1" );
}

sub login_cas {
	my $self        = shift;
	my %opts        = @_;
	my $cookie_name = $opts{cookie_name};
	my $ticket      = $opts{ticket};        #cas ticket from iPlant
	my $this_url    = $opts{this_url};      #not sure what this does
	my $coge        = $opts{coge};          #coge object

	#print STDERR Dumper \%opts;
	my $ua = new LWP::UserAgent;

	my $request = '<SOAP-ENV:Envelope xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/">'.
		'<SOAP-ENV:Header/><SOAP-ENV:Body><samlp:Request xmlns:samlp="urn:oasis:names:tc:SAML:1.0:protocol"  MajorVersion="1" MinorVersion="1" RequestID="_192.168.167.84.1024506224022"  IssueInstant="2010-05-13T16:43:48.099Z"><samlp:AssertionArtifact>'
	  . $ticket . '</samlp:AssertionArtifact></samlp:Request></SOAP-ENV:Body></SOAP-ENV:Envelope>';

	my $request_ua = HTTP::Request->new( POST => 'https://auth.iplantcollaborative.org/cas/samlValidate?TARGET=' . $this_url );
	$request_ua->content($request);
	$request_ua->content_type("text/xml; charset=utf-8");
	my $response = $ua->request($request_ua);
	my $result   = $response->content;
	my $uname;
	my $fname;
	my $lname;
	my $email;

	if ($result) {
		( $uname, $fname, $lname, $email ) = parse_saml_response($result);
	}
	return unless $uname; # Not logged in.  Return
	
	my $coge_user;
	($coge_user) = $coge->resultset('User')->search( { user_name => $uname } );
	unless ($coge_user) {
		# Create new user
		$coge_user = $coge->resultset('User')->create(
			{	user_name   => $uname,
				first_name  => $fname,
				last_name   => $lname,
				email       => $email,
				description => "Validated by iPlant"
			}
		); #do we have a valid user in the database, if not create
		$coge_user->insert;

		$coge->resultset('Log')->create( { user_id => $coge_user->id, page => 'Web.pm', description => 'create user' } );
	}

	#create a session ID for the user and log
	my $session = md5_base64( $uname . $ENV{REMOTE_ADDR} );
	$session =~ s/\+/1/g;
	my $sid = $coge->log_user( user => $coge_user, session => $session );

	# mdb added 10/19/12 - FIXME key/secret are hardcoded - wait: this will get replaced by openauth soon
	#$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME}=0; # this doesn't work for bypassing cert check, need line in apache cfg
	$request_ua = HTTP::Request->new( POST => 'https://user.iplantcollaborative.org/api/v1/service/coge/add/' . $uname );
	$request_ua->authorization_basic('6mv9x9lyts8oje8uj3t6yo', 'f59ba33ee35d363ffefd8b27b375e587b0e5c7a1');
	#print STDERR "request uri: " . $request_ua->uri . "\n";
	#$request_ua->content($request);
	$request_ua->content_type("text/xml; charset=utf-8");
	$response = $ua->request($request_ua);
	#if ($response->is_success()) {
		#print STDERR "status_line: " . $response->status_line() . "\n";
		#my $header = $response->header;
		$result = $response->content;
		#print STDERR "content: <begin>$result<end>\n";
	#}
	#else {
	#	print STDERR "bad response\n";	
	#}

	#gen and set the web cookie, yum!
	my $c = CoGe::Accessory::LogUser->gen_cookie(
		session     => $session,
		cookie_name => $cookie_name,
	);

#	print STDERR "login_cas:  gen_cookie " . (Dumper $c) . "\n";
	print CGI::header( -cookie => [$c] );
	return $coge_user;
}

sub parse_saml_response {
	my $response = $_[0];

	# mdb modified 4/4/13 for iPlant CAS update - XML::Simple doesn't support namespaces
	if ( $response =~ m/saml1p:Success/ ) {
		my $ref = XMLin($response);
		my ($user_id) =
		  $ref->{'SOAP-ENV:Body'}->{'saml1p:Response'}->{'saml1:Assertion'}->{'saml1:AttributeStatement'}
		  ->{'saml1:Subject'}->{'saml1:NameIdentifier'};
		my @tmp = @{$ref->{'SOAP-ENV:Body'}->{'saml1p:Response'}->{'saml1:Assertion'}->{'saml1:AttributeStatement'}->{'saml1:Attribute'}};
		my %attr =
		  map { $_->{'AttributeName'}, $_->{'saml1:AttributeValue'} }
		  @{ $ref->{'SOAP-ENV:Body'}->{'saml1p:Response'}->{'saml1:Assertion'}->{'saml1:AttributeStatement'}->{'saml1:Attribute'} };
		my ($user_lname) = $attr{lastName}->{content};
		my ($user_fname) = $attr{firstName}->{content};
		my ($user_email) = $attr{email}->{content};

		#print STDERR "parse_saml_response: ".$user_id.'   '.$user_fname.'   '.$user_lname.'  '.$user_email."\n";
		return ( $user_id, $user_fname, $user_lname, $user_email );
	}
}

sub login { # FIXME mdb 9/21/12 - what is this?
	my %opts  = @_;
	my $coge  = $opts{coge};
	my $uname = $opts{uname};
	my $fname = $opts{fname};
	my $lname = $opts{lname};
	my $email = $opts{email};
	my $url   = $opts{url};

}

sub ajax_func {
	return (
		read_log            => \&read_log,
		initialize_basefile => \&initialize_basefile,
	);
}

sub log_history {
	my %opts = @_;
	my $db				= $opts{db};
	my $user_id 		= $opts{user_id};
	my $page 			= $opts{page};
	my $description 	= $opts{description};
	my $link			= $opts{link};
	
	$user_id = 0 unless (defined $user_id);
	$db->resultset('Log')->create( { user_id => $user_id, page => $page, description => $description, link => $link } );
}

sub get_tiny_link {
	my %opts = @_;
	my $url 			= $opts{url};
	my $db				= $opts{db};
	my $user_id 		= $opts{user_id};
	my $page 			= $opts{page};
	my $log_msg			= $opts{log_msg};
	my $disable_logging	= $opts{disable_logging}; # flag
	$url =~ s/:::/__/g;
	my $tiny = LWP::Simple::get("http://genomevolution.org/r/yourls-api.php?signature=d57f67d3d9&action=shorturl&format=simple&url=$url");
	unless ($tiny) {
		return "Unable to produce tiny url from server";
	}

	# Log the page
	if ($db and not $disable_logging) {
		$page =~ s/.pl$//; # remove trailing .pl extension
		log_history( db => $db, user_id => $user_id, page => $page, description => ($log_msg ? $log_msg : 'page access'), link => $tiny );
	}

	return $tiny;
}

sub get_job {
    my %args = @_;
    my $job;
    my $tiny_link = $args{tiny_link};
    my $user_id = $args{user_id};
    my $title = $args{title};
    my $coge = $args{db_object};


    $user_id = 0 unless defined($user_id);

    my $prev_submission = $coge->resultset('Job')->search({
        user_id => $user_id,
        link => $tiny_link
    });

    if($prev_submission->count < 1) {
            $job = $coge->resultset('Job')->create({
                "link" => $tiny_link,
                "page" => $title,
                "process_id" => getpid(),
                "user_id" => $user_id,
                "status" => 1
            });
    } else {
        $job = $prev_submission->next;
        $job->update({
            status => 1,
            process_id => getpid()
        });
    }

    return $job;
}

sub write_log {
	$| = 1;
	my $message = shift;
	$message =~ /(.*)/xs;
	$message = $1;
	my $file = shift;
	return unless $file;
	open( OUT, ">>$file" ) || return;
	print OUT $message, "\n";
	close OUT;
}

sub read_log {
	my %args    = @_;
	my $logfile = $args{logfile};
	my $prog    = $args{prog};
	my $tempdir = $args{tempdir};
	$tempdir = $TEMPDIR unless $tempdir;
	return unless $logfile;
	$logfile .= ".log" unless $logfile =~ /log$/;
	unless ( $logfile =~ /^$tempdir/ ) {
		$logfile = "$prog/" . $logfile if $prog;
		$logfile = "$tempdir/" . $logfile;
	}
	return unless -r $logfile;
	my $str;
	open( IN, $logfile );
	while (<IN>) {
		$str .= $_;
	}
	close IN;
	return $str;
}

sub check_filename_taint {
	my $v = shift;
	return 1 unless $v;
	if ( $v =~ /^([A-Za-z0-9\-\.=\/_#]*)$/ ) {
		my $v1 = $1;
		return ($v1);
	}
	else {
		return (0);
	}
}

sub check_taint {
	my $v = shift;
	return 1 unless $v;
	if ( $v =~ /^([-\w\._=\s+\/,#\]\['"%]+)$/ ) {
		$v = $1;

		# $v now untainted
		return ( 1, $v );
	}
	else {

		# data should be thrown out
		carp "'$v' failed taint check\n";
		return (0);
	}
}

sub save_settings {
	my %opts    = @_;
	my $user    = $opts{user};
	my $user_id = $opts{user_id};
	my $page    = $opts{page};
	my $opts    = $opts{opts};
	my $coge    = $opts{coge};
	$opts = Dumper $opts unless $opts =~ /VAR1/;
	$user_id = $user->id if ( ref($user) =~ /User/i ) && !$user_id;

	unless ($user_id) {
		my ($user_obj) = $coge->resultset('User')->search( { user_name => $user } );
		$user_id = $user_obj->id if $user_obj;
	}
	return unless $user_id;

	#delete previous settings
	foreach my $item ( $coge->resultset('WebPreferences')->search( { user_id => $user_id, page => $page } ) )
	{
		$item->delete;
	}
	my $item = $coge->resultset('WebPreferences')->new( { user_id => $user_id, page => $page, options => $opts } );
	$item->insert;
	return $item;
}

sub load_settings {
	my %opts    = @_;
	my $user    = $opts{user};
	my $user_id = $opts{user_id};
	my $page    = $opts{page};
	my $coge    = $opts{coge};
	unless ($coge) {
		print STDERR "need a valid coge object";
		return;
	}
	$user_id = $user->id if ( ref($user) =~ /User/i ) && !$user_id;
	unless ($user_id) {
		my ($user_obj) = $coge->resultset('User')->search( { user_name => $user } );
		$user_id = $user_obj->id if $user_obj;
	}
	return {} unless $user_id;
	my ($item) = $coge->resultset('WebPreferences')->search( { user_id => $user_id, page => $page } );
	return {} unless $item;
	my $prefs;
	my $opts = $item->options if $item;
	return {} unless $opts;
	$opts =~ s/VAR1/prefs/;
	eval $opts;
	return $prefs;
}

sub reset_settings {
	my %opts    = @_;
	my $user    = $opts{user};
	my $user_id = $opts{user_id};
	my $page    = $opts{page};
	my $coge    = $opts{coge};
	$user_id = $user->id if ( ref($user) =~ /User/i ) && !$user_id;
	unless ($user_id) {
		my ($user_obj) = $coge->resultset('User')->search( { user_name => $user } );
		$user_id = $user_obj->id if $user_obj;
	}
	return unless $user_id;
	my ($item) = $coge->resultset('WebPreferences')->search( { user_id => $user_id, page => $page } );
	$item->delete;
}

sub initialize_basefile {
	my ( $self, %opts ) = self_or_default(@_);
	my $basename    = $opts{basename};
	my $prog        = $opts{prog};
	my $return_name = $opts{return_name};
	my $tempdir     = $opts{tempdir} || $TEMPDIR;
	$tempdir .= "/" . $prog if $prog;
	if ($basename) {
		#print STDERR "Have basename: $basename\n";
		($basename) = $basename =~ /([^\/].*$)/;
		my ( $x, $cleanname ) = check_taint($basename);
		$self->basefilename($cleanname);
		my $basefile = $tempdir . "/" . $cleanname;
		$basefile =~ s/\/\/+/\//g;
		$self->basefile($basefile);
		$self->logfile( $self->basefile . ".log" );
		$self->sqlitefile( $self->basefile . ".sqlite" );
	}
	else {
		mkdir "$tempdir", 0777 unless -d "$tempdir";
		$prog = "CoGe" unless $prog;
		my $file = new File::Temp(
			TEMPLATE => $prog . '_XXXXXXXX',
			DIR      => "$tempdir/",
			#SUFFIX=>'.png',
			UNLINK => 1
		);
		$self->basefile( $file->filename );
		$self->logfile( $self->basefile . ".log" );
		$self->sqlitefile( $self->basefile . ".sqlite" );
		$self->basefilename( $file->filename =~ /([^\/]*$)/ );
	}

#    print STDERR "Basename: ",$self->basefilename,"\n";
#    print STDERR "sqlitefile: ",$self->sqlitefile,"\n";
#    print STDERR "Basefile: ",$self->basefile,"\n";

	if ( -r $self->logfile && !$basename ) {
		print STDERR "in Web.pm sub initialize_basefile.  Logfile " . $self->logfile . " already exist.  Possible problem.  Regenerating basefile.\n";
		return $self->initialize_basefile(%opts);
	}
	elsif ($return_name) {
		return $self->basefilename;
	}
	else { return $self; }
}

sub gzip {
	my ( $self, $file, $conf_file ) = self_or_default(@_);
	$conf_file = $ENV{HOME} . 'coge.conf' unless $conf_file;
	my $P    = $self->get_defaults($conf_file);
	my $GZIP = $P->{GZIP};
	return $file unless $file;
	return $file . ".gz" if -r "$file.gz";
	return $file unless -r $file;
	return $file if $file =~ /\.gz$/;
	`$GZIP $file` if -r $file;
	my $tmp = $file . ".gz";
	return -r $tmp ? $tmp : $file;
}

sub gunzip {
	my ( $self, $file, $conf_file, $debug ) = self_or_default(@_);
	$conf_file = $ENV{HOME} . 'coge.conf' unless $conf_file;
	my $P      = $self->get_defaults($conf_file);
	my $GUNZIP = $P->{GUNZIP};
	unless ($GUNZIP) {
		print STDERR "ERROR: in gunzip!  gunzip binary is not specified!\n"
		  if $debug;
	}
	print STDERR "Debugging sub gunzip\n" if $debug;
	print STDERR "\t", $file, "\n" if $debug;
	$file .= ".gz" if -r $file . ".gz";
	print STDERR "\t", $file, "!\n" if $debug;
	if ( -r $file && $file =~ /\.gz/ ) {
		print STDERR "\t", "Running $GUNZIP $file\n";
		`$GUNZIP $file`;
	}
	my $tmp = $file;
	$tmp =~ s/\.gz$//;
	my $return = -r $tmp ? $tmp : $file;
	print STDERR "\t", "returning $return\n" if $debug;
	return $return;
}

sub irods_ils {
	my $path = shift;
	$path = '' unless $path;
#	print STDERR "irods_ils: path=$path\n";

	my $P = get_defaults( $ENV{HOME} . 'coge.conf' );
	my $env_file = $P->{IRODSENV};
	if (not defined $env_file or not -e $env_file) {
		print STDERR "fatal error: iRODS env file missing!\n";
		return { error => "Error: iRODS env file missing" };
	}

	my $cmd = "export irodsEnvFile='$env_file'; ils -l $path 2>&1";
#	print STDERR "cmd: $cmd\n";
#	my @ils = `$cmd`; # old way of executing command, replaced by better error checking below
	my @ils = capture(EXIT_ANY, $cmd);
   	if ($EXITVAL) {
		return { error => "Error: ils rc=$EXITVAL" };
	}

	$path = shift @ils;
#	if ($path =~ /^ERROR/) { # iRODS error message
#		my $result = { type => 'error', name => $path };
#		return wantarray ? ($result) : [$result];
#	}
	chomp($path);
	chop($path);

	my @result;
	foreach my $line (@ils) {
		my ($type, $size, $timestamp, $name);

		chomp $line;
		if ($line =~ /^\s*C\-/) { # directory
			$type = 'directory';
			($name) = $line =~ /([^\/\s]+)\s*$/;
			if ($name) { $name .= '/'; }
			else { $name = 'error' };
			($size, $timestamp) = ('', '');
		}
		else { # file
			$type = 'file';
			(undef, undef, undef, undef, $size, $timestamp, undef, $name) = split(/\s+/, $line);
		}
		
		push @result, 
			{ type => $type, 
			  size => $size, 
			  timestamp => $timestamp, 
			  name => $name,
			  path => $path . '/' . $name
			};
	}
	@result = sort {$a->{type} cmp $b->{type}} @result; # directories before files
	
	return { items => \@result };
}

sub irods_chksum {
	my $path = shift;
	return 0 unless ($path);

	my $P = get_defaults( $ENV{HOME} . 'coge.conf' );
	my $env_file = $P->{IRODSENV};
	if (not defined $env_file or not -e $env_file) {
		print STDERR "fatal error: iRODS env file missing!\n";
		return;	
	}

	my $cmd = "export irodsEnvFile='$env_file'; ichksum $path";
#	print STDERR "cmd: $cmd\n";
	my @output = `$cmd`;
	my ($chksum) = $output[0] =~ /\s*\S+\s+(\S+)/;
#	print STDERR "chksum: $chksum\n";

	return $chksum;
}

sub irods_iget {
	my ($src, $dest) = @_;
#	print STDERR "irods_iget $src $dest\n";

	my $P = get_defaults( $ENV{HOME} . 'coge.conf' );
	my $env_file = $P->{IRODSENV};
	if (not defined $env_file or not -e $env_file) {
		print STDERR "fatal error: iRODS env file missing!\n";
		return;	
	}

	my $cmd = "export irodsEnvFile='$env_file'; iget -fT $src $dest";
#	print STDERR "cmd: $cmd\n";
	my @ils = `$cmd`;
#	print STDERR "@ils";
	
	return;
}

sub send_email {
	my %opts = @_;
	my $from 	= $opts{from};
	my $to 		= $opts{to};
	my $subject = $opts{subject};
	my $body 	= $opts{body};

	print STDERR "Sending email: $from $to $subject\n";

	my $mailer = Mail::Mailer->new("sendmail");
	$mailer->open({
		From    => $from,
		To      => $to,
		Subject => $subject,
	}) or die "Can't open: $!\n";

	print $mailer $body;
	$mailer->close();
}

1;


=head1 NAME

Web

=head1 SYNOPSIS

use Web

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
