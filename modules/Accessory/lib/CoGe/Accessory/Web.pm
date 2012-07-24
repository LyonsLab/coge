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

BEGIN {
  use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK $Q $cogex $TEMPDIR $BASEDIR);
    require Exporter;

    $BASEDIR="/opt/apache/CoGe/";
    $VERSION     = 0.1;
    $TEMPDIR = $BASEDIR."tmp";
    @ISA         = (@ISA, qw (Exporter));

    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();#qw (login write_log read_log check_taint check_filename_taint save_settings load_settings reset_settings initialize_basefile);
#    $cogex = CoGeX->dbconnect();
#    $cogex->storage->debugobj(new DBIxProfiler());
#    $cogex->storage->debug(1);
    __PACKAGE__->mk_accessors qw(restricted_orgs basefilename basefile logfile sqlitefile); 
  
}

sub get_defaults
  {
    my ($self, $param_file) = self_or_default(@_);
    $param_file = $BASEDIR."coge.conf" unless defined $param_file;
    unless (-r $param_file)
      {
	print STDERR qq{Either no parameter file specified or unable to read paramer file ($param_file).
A valid parameter file must be specified or very little will work!};
	return 0;
      }
    open (IN, $param_file);
    my %items;
    while (<IN>)
      {
	chomp;
	next if /^#/;
	next unless $_;
	my ($name, $path) = split(/\s+/,$_,2);
	$items{$name}=$path;
      }
    close IN;
    return \%items;
  }

sub dataset_search_for_feat_name
  {
    my ($self, $accn, $num, $dsid, $featid, $coge) = self_or_default(@_);
    $num = 1 unless $num;
    return ( qq{<input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">}, $num )unless $accn;
    my $html;
    my %sources;
    my %restricted_orgs = %{$self->restricted_orgs} if $self->restricted_orgs;
    my $rs = $coge->resultset('Dataset')->search(
						  {
						   'feature_names.name'=> $accn,
						  },
						  {
						   'join'=>{
							    'features' => 'feature_names',
							   },
							    
						   'prefetch'=>['datasource', 'organism'],
						  }
						 );
    while (my $ds = $rs->next())
      {
	my $name = $ds->name;
	my $ver = $ds->version;
	my $desc = $ds->description;
	my $sname = $ds->datasource->name;
	my $ds_name = $ds->name;
	my $org = $ds->organism->name;
	my $title = "$org: $ds_name ($sname, v$ver)";
	next if $restricted_orgs{$org};
	$sources{$ds->id} = {
			     title=>$title,
			     version=>$ver,
			    };
      }
     if (keys %sources)
       {
 	$html .= qq{
 <SELECT name = "dsid$num" id= "dsid$num" onChange="feat_search(['accn$num','dsid$num', 'args__$num'],['feat$num']);" >
 };
 	foreach my $id (sort {$sources{$b}{version} <=> $sources{$a}{version}} keys %sources)
 	  {
 	    my $val = $sources{$id}{title};
 	    $html  .= qq{  <option value="$id"};
	    $html .= qq{ selected } if $dsid && $id == $dsid;
	    $html .= qq{>$val\n};
 	  }
 	$html .= qq{</SELECT>\n};
 	my $count = scalar keys %sources;
 	$html .= qq{<font class=small>($count)</font>};
       }
     else
       {
 	$html .= qq{Accession not found <input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">\n};	
       }    
    return ($html,$num);
  }

sub feat_search_for_feat_name
  {
    my ($self, $accn, $dsid, $num, $coge) = self_or_default(@_);
    return qq{<input type="hidden" id="featid$num">\n} unless $dsid;
    my @feats;
    my $rs = $coge->resultset('Feature')->search(
						  {
						   'feature_names.name'=> $accn,
						   'dataset.dataset_id' => "$dsid",
						  },
						  {
						   'join'=>['feature_type','dataset', 'feature_names'],
						   'prefetch'=>['feature_type', 'dataset'],
						  }
						 );
    my %seen;
    while( my $f =$rs->next())
      {
	next unless $f->dataset->id == $dsid;
#	next if $f->feature_type->name =~ /CDS/i;
#	next if $f->feature_type->name =~ /RNA/i;
	push @feats, $f unless $seen{$f->id};
	$seen{$f->id}=1;
      }
    my $html;
    if (@feats)
      {
	$html .= qq{
<SELECT name = "featid$num" id = "featid$num" >
  };
	foreach my $feat (sort {$a->type->name cmp $b->type->name} @feats)
	  {
	    my $loc = "(".$feat->type->name.") Chr:".$feat->locations->next->chromosome." ".$feat->start."-".$feat->stop;
	    #working here, need to implement genbank_location_string before I can progress.  Need 
	    $loc =~ s/(complement)|(join)//g;
	    my $fid = $feat->id;
	    $html .= qq {  <option value="$fid">$loc \n};
	  }
	$html .= qq{</SELECT>\n};
	my $count = scalar @feats;
	$html .= qq{<font class=small>($count)</font>};
      }
    else
      {
	$html .=  qq{<input type="hidden" id="featid$num">\n}
      }
    return $html;
  }

sub self_or_default { #from CGI.pm
    return @_ if defined($_[0]) && (!ref($_[0])) &&($_[0] eq 'CoGe::Accessory::Web');
    unless (defined($_[0]) && 
            (ref($_[0]) eq 'CoGe::Accessory::Web' || UNIVERSAL::isa($_[0],'CoGe::Accessory::Web')) # slightly optimized for common case
            ) {
        $Q = CoGe::Accessory::Web->new unless defined($Q);
        unshift(@_,$Q);
    }
    return wantarray ? @_ : $Q;
}

sub logout_cas {
  my $self = shift;
  my %opts = @_;
  my $cookie_name = $opts{cookie_name};
  my $coge = $opts{coge};
  my $user = $opts{user};
  my $form = $opts{form}; #CGI form for calling page
  my $url = $opts{this_url};
  my $cas_url ='http://coge.iplantcollaborative.org/coge/';
  my %cookies = fetch CGI::Cookie;
  $url = $form->url() unless $url;
  my $session = md5_base64($user->user_name.$ENV{REMOTE_ADDR});
  $session =~ s/\+/1/g;
  ($session) = $coge->resultset('UserSession')->find({session=>$session});
  $session->delete if $session;
  print "Location: ".$form->redirect("https://auth.iplantcollaborative.org/cas/logout?service=".$url."&gateway=1");
}

sub login_cas{
  my $self = shift;
  my %opts = @_;
  my $cookie_name = $opts{cookie_name};
  my $ticket = $opts{ticket}; #cas ticket from iPlant
  my $this_url = $opts{this_url}; #not sure what this does
  my $coge = $opts{coge}; #coge object
#  print STDERR Dumper \%opts;
  my $ua = new LWP::UserAgent;

  my $request = '<SOAP-ENV:Envelope xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"><SOAP-ENV:Header/><SOAP-ENV:Body><samlp:Request xmlns:samlp="urn:oasis:names:tc:SAML:1.0:protocol"  MajorVersion="1" MinorVersion="1" RequestID="_192.168.167.84.1024506224022"  IssueInstant="2010-05-13T16:43:48.099Z"><samlp:AssertionArtifact>'.$ticket.'</samlp:AssertionArtifact></samlp:Request></SOAP-ENV:Body></SOAP-ENV:Envelope>';
  
  my $request_ua = HTTP::Request->new(POST => 'https://auth.iplantcollaborative.org/cas/samlValidate?TARGET='.$this_url);
  $request_ua->content($request);
  $request_ua->content_type("text/xml; charset=utf-8");
  my $response = $ua->request($request_ua);
  my $result =$response->content;
  my $uname;
  my $fname;
  my $lname;
  my $email;
  if($result){
    ($uname,$fname,$lname,$email) = parse_saml_response($result);
  }
  return unless $uname; #not logged in.  Return
  my $coge_user;
  ($coge_user) = $coge->resultset('User')->search({user_name=>$uname});
  unless ($coge_user)
    {
      $coge_user =  $coge->resultset('User')->create({user_name=>$uname,first_name=>$fname,last_name=>$lname,email=>$email, description=>"Validated by iPlant"}); #do we have a valid user in the database, if not create
      $coge_user->insert;
    }
  #create a session ID for the user and log
  my $session = md5_base64($uname.$ENV{REMOTE_ADDR});
  $session =~ s/\+/1/g;
  my $sid = $coge->log_user(user=>$coge_user,session=>$session);
  #gen and set the web cookie, yum!
  my $c = CoGe::Accessory::LogUser->gen_cookie(session=>$session,cookie_name=>$cookie_name);
#  print "Set-Cookie: $c\n";
  print CGI::header(-cookie=>[$c]);
  return $coge_user;
}

sub parse_saml_response{
	
	my $response = $_[0];
	
	if($response =~ m/samlp:Success/){
		
		my $ref = XMLin($response);
		my ($user_id) = $ref->{'SOAP-ENV:Body'}->{Response}->{Assertion}->{AttributeStatement}->{Subject}->{NameIdentifier};
		my %attr = map {$_->{'AttributeName'}, $_->{'AttributeValue'}} @{$ref->{'SOAP-ENV:Body'}->{Response}->{Assertion}->{AttributeStatement}->{Attribute}};
		my ($user_lname) = $attr{lastName};
		my ($user_fname) = $attr{firstName};
		my ($user_email) = $attr{email};
#		print STDERR $user_id.'   '.$user_fname.'   '.$user_lname.'  '.$user_email;
		
		return ($user_id,$user_fname,$user_lname,$user_email);
	}
}

sub login 
  {
    my %opts = @_;
    my $coge = $opts{coge};
    my $uname = $opts{uname};
    my $fname = $opts{fname};
    my $lname = $opts{lname};
    my $email = $opts{email};
    my $url = $opts{url};
    
  }
	   
sub ajax_func
  {
    return 
      (
       read_log=>\&read_log,
       initialize_basefile=>\&initialize_basefile,
      );
  }


sub get_tiny_link
  {
    my %opts = @_;
    my $url = $opts{url};
    $url =~ s/:::/__/g;
    my $html;
    my $tiny = LWP::Simple::get("http://genomevolution.org/r/yourls-api.php?signature=d57f67d3d9&action=shorturl&format=simple&url=$url");
    unless ($tiny)
      {
        return "Unable to produce tiny url from server";
      }
    return $tiny;
  }



sub write_log
  {
    $| = 1;
    my $message = shift;
    $message =~ /(.*)/xs;
    $message = $1;
    my $file = shift;
    return unless $file;
    open (OUT, ">>$file") || return;
    print OUT $message,"\n";
    close OUT;
  }

sub read_log
  {
    my %args = @_;
    my $logfile = $args{logfile};
    my $prog = $args{prog};
    my $tempdir = $args{tempdir};
    $tempdir = $TEMPDIR unless $tempdir;
    return unless $logfile;
    $logfile .= ".log" unless $logfile =~ /log$/;
    unless ($logfile =~ /^$tempdir/)
      {
	$logfile = "$prog/".$logfile if $prog;
	$logfile = "$tempdir/".$logfile;
      }
    return unless -r $logfile;
    my $str;
    open (IN, $logfile);
    while (<IN>)
      {
	$str .= $_;
      }
    close IN;
    return $str;
  }

sub check_filename_taint {
  my $v = shift;
  return 1 unless $v;
  if ($v =~ /^([A-Za-z0-9\-\.=\/_#]*)$/) {
    my $v1 = $1;
    return($v1);
  } else {
    return(0);
  }
}

sub check_taint {
  my $v = shift;
  return 1 unless $v;
  if ($v =~ /^([-\w\._=\s+\/,#\]\['"%]+)$/) {
    $v = $1;
    # $v now untainted
    return(1,$v);
  } else {
    # data should be thrown out
    carp "'$v' failed taint check\n";
    return(0);
  }
}

sub save_settings
  {
    my %opts = @_;
    my $user = $opts{user};
    my $user_id = $opts{user_id};
    my $page = $opts{page};
    my $opts = $opts{opts};
    my $coge = $opts{coge};
    $opts = Dumper $opts unless $opts =~ /VAR1/;
    $user_id = $user->id if (ref ($user) =~ /User/i) && !$user_id;
    unless ($user_id)
      {
	my ($user_obj) = $coge->resultset('User')->search({user_name=>$user});
	$user_id = $user_obj->id if $user_obj;
      }
    return unless $user_id;
    #delete previous settings
    foreach my $item ($coge->resultset('WebPreferences')->search({user_id=>$user_id, page=>$page}))
      {
	$item->delete;
      }
    my $item = $coge->resultset('WebPreferences')->new({user_id=>$user_id, page=>$page, options=>$opts});
    $item->insert;
    return $item;
  }

sub load_settings
  {
    my %opts = @_;
    my $user = $opts{user};
    my $user_id = $opts{user_id};
    my $page = $opts{page};
    my $coge = $opts{coge};
    unless ($coge)
     {
 	print STDERR "need a valid coge object";
	return;
     }
    $user_id = $user->id if (ref ($user) =~ /User/i) && !$user_id;
    unless ($user_id)
      {
	my ($user_obj) = $coge->resultset('User')->search({user_name=>$user});
	$user_id = $user_obj->id if $user_obj;
      }
    return {} unless $user_id;
    my ($item) = $coge->resultset('WebPreferences')->search({user_id=>$user_id, page=>$page});
    return {} unless $item;
    my $prefs;
    my $opts = $item->options if $item;
    return {} unless $opts;
    $opts =~ s/VAR1/prefs/;
    eval $opts;
    return $prefs;
  }

sub reset_settings
  {
    my %opts = @_;
    my $user = $opts{user};
    my $user_id = $opts{user_id};
    my $page = $opts{page};
    my $coge = $opts{coge};
    $user_id = $user->id if (ref ($user) =~ /User/i) && !$user_id;
    unless ($user_id)
      {
	my ($user_obj) = $coge->resultset('User')->search({user_name=>$user});
	$user_id = $user_obj->id if $user_obj;
      }
    return unless $user_id;
    my ($item) = $coge->resultset('WebPreferences')->search({user_id=>$user_id, page=>$page});
    $item->delete;
  }

sub initialize_basefile
  {

    my ($self, %opts) = self_or_default(@_);
    my $basename = $opts{basename};
    my $prog=$opts{prog};
    my $return_name = $opts{return_name};
    my $tempdir = $opts{tempdir} || $TEMPDIR;
    $tempdir .= "/".$prog if $prog;
    if ($basename)
      {
#	print STDERR "Have basename: $basename\n";
	($basename) = $basename =~ /([^\/].*$)/;
	my ($x, $cleanname) = check_taint($basename);
	$self->basefilename($cleanname);
	my $basefile = $tempdir."/".$cleanname;
	$basefile =~ s/\/\/+/\//g;
	$self->basefile($basefile);
	$self->logfile($self->basefile.".log");
	$self->sqlitefile($self->basefile.".sqlite");
      }
    else
      {
	mkdir "$tempdir",0777 unless -d "$tempdir";
	$prog = "CoGe" unless $prog;
	my $file = new File::Temp ( TEMPLATE=>$prog.'_XXXXXXXX',
				    DIR=>"$tempdir/",
				    #SUFFIX=>'.png',
				    UNLINK=>1);
	$self->basefile($file->filename);
	$self->logfile($self->basefile.".log");
	$self->sqlitefile($self->basefile.".sqlite");
	$self->basefilename($file->filename =~ /([^\/]*$)/)
      }
#    print STDERR "Basename: ",$self->basefilename,"\n";
#    print STDERR "sqlitefile: ",$self->sqlitefile,"\n";
#    print STDERR "Basefile: ",$self->basefile,"\n";
    
    if (-r $self->logfile && ! $basename)
      {
	print STDERR "in Web.pm sub initialize_basefile.  Logfile ".$self->logfile." already exist.  Possible problem.  Regenerating basefile.\n";
	return $self->initialize_basefile(%opts);
      }
    elsif ($return_name)
      {
	return $self->basefilename;
      }
    else {return $self;}
  }


sub gzip
    {
      my ($self, $file, $conf_file) = self_or_default(@_);
      $conf_file = $ENV{HOME}.'coge.conf' unless $conf_file;
      my $P = $self->get_defaults($conf_file);
      my $GZIP = $P->{GZIP};
      return $file unless $file;
      return $file.".gz" if -r "$file.gz";
      return $file unless -r $file;
      return $file if $file =~ /\.gz$/;
      `$GZIP $file` if -r $file;
      my $tmp = $file.".gz";
      return -r $tmp ? $tmp : $file;
    }

sub gunzip
    {
      my ($self, $file, $conf_file, $debug) = self_or_default(@_);
      $conf_file = $ENV{HOME}.'coge.conf' unless $conf_file;
      my $P = $self->get_defaults($conf_file);
      my $GUNZIP = $P->{GUNZIP};
      unless ($GUNZIP)
	{
	  print STDERR "ERROR: in gunzip!  gunzip binary is not specified!\n" if $debug;
	}
      print STDERR "Debugging sub gunzip\n" if $debug;
      print STDERR "\t",$file,"\n" if $debug;
      $file .= ".gz" if -r $file .".gz";
      print STDERR "\t",$file,"!\n" if $debug;
      if (-r $file && $file =~ /\.gz/)
	{
	  print STDERR "\t", "Running $GUNZIP $file\n";
	  `$GUNZIP $file`;
	}
      my $tmp = $file;
      $tmp =~ s/\.gz$//;
      my $return = -r $tmp ? $tmp : $file;
      print STDERR "\t","returning $return\n" if $debug;
      return $return;
    }


1;
