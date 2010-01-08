package CoGe::Accessory::Web;

use strict;
use CoGeX;
use Data::Dumper;
use base 'Class::Accessor';
use CGI::Carp('fatalsToBrowser');
use CGI;
use DBIxProfiler;
use File::Basename;
use File::Temp;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $Q $cogex $TEMPDIR);
    $VERSION     = 0.1;
    $TEMPDIR = "/opt/apache/CoGe/tmp";
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw (login write_log read_log check_taint check_filename_taint save_settings load_settings reset_settings initialize_basefile);
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    $cogex = CoGeX->dbconnect();
#    $cogex->storage->debugobj(new DBIxProfiler());
#    $cogex->storage->debug(1);
    __PACKAGE__->mk_accessors qw(restricted_orgs basefilename basefile logfile sqlitefile);
 }

sub dataset_search_for_feat_name
  {
    my ($self, $accn, $num, $dsid, $featid) = self_or_default(@_);
    $num = 1 unless $num;
    return ( qq{<input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">}, $num )unless $accn;
    my $html;
    my %sources;
    my %restricted_orgs = %{$self->restricted_orgs} if $self->restricted_orgs;
    my $rs = $cogex->resultset('Dataset')->search(
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
    my ($self, $accn, $dsid, $num) = self_or_default(@_);
    return qq{<input type="hidden" id="featid$num">\n} unless $dsid;
    my @feats;
    my $rs = $cogex->resultset('Feature')->search(
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

sub login
  {
    my $form = new CGI;
    my $url = "index.pl?url=".$form->url(-relative=>1, -query=>1);
#    print STDERR $url;
    $url =~ s/&|;/:::/g;
    my $html1 = qq{
<SCRIPT language="JavaScript">
window.location=$url;
</SCRIPT>
};
    my $html = qq{
<html>
<head>
<title>CoGe:  the best program of its kind, ever.</title>
<meta http-equiv="REFRESH" content="1;url=$url"></HEAD>
<BODY>
You are not logged in.  You will be redirected to <a href = $url>here</a> in one second.
</BODY>
</HTML>
};
    return $html;
  }

sub ajax_func
  {
    return 
      (
       login=>\&login,
       read_log=>\&read_log,
       initialize_basefile=>\&initialize_basefile,
      );
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
#    print STDERR "Checking logfile $logfile\n";
    return unless $logfile;
    $logfile .= ".log" unless $logfile =~ /log$/;
    $logfile = $TEMPDIR."/$prog/$logfile" unless $logfile =~ /^$TEMPDIR/;
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
  if ($v =~ /^([A-Za-z0-9\-\.=\/_]*)$/) {
    my $v1 = $1;
    return($v1);
  } else {
    return(0);
  }
}

sub check_taint {
  my $v = shift;
  return 1 unless $v;
  if ($v =~ /^([-\w._=\s+\/]+)$/) {
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
    $opts = Dumper $opts unless $opts =~ /VAR1/;
    $user_id = $user->id if (ref ($user) =~ /User/i) && !$user_id;
    unless ($user_id)
      {
	my ($user_obj) = $cogex->resultset('User')->search({user_name=>$user});
	$user_id = $user_obj->id if $user_obj;
      }
    return unless $user_id;
    #delete previous settings
    foreach my $item ($cogex->resultset('WebPreferences')->search({user_id=>$user_id, page=>$page}))
      {
	$item->delete;
      }
    my $item = $cogex->resultset('WebPreferences')->new({user_id=>$user_id, page=>$page, options=>$opts});
    $item->insert;
    return $item;
  }

sub load_settings
  {
    my %opts = @_;
    my $user = $opts{user};
    my $user_id = $opts{user_id};
    my $page = $opts{page};
    $user_id = $user->id if (ref ($user) =~ /User/i) && !$user_id;
    unless ($user_id)
      {
	my ($user_obj) = $cogex->resultset('User')->search({user_name=>$user});
	$user_id = $user_obj->id if $user_obj;
      }
    return {} unless $user_id;
    my ($item) = $cogex->resultset('WebPreferences')->search({user_id=>$user_id, page=>$page});
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
    $user_id = $user->id if (ref ($user) =~ /User/i) && !$user_id;
    unless ($user_id)
      {
	my ($user_obj) = $cogex->resultset('User')->search({user_name=>$user});
	$user_id = $user_obj->id if $user_obj;
      }
    return unless $user_id;
    my ($item) = $cogex->resultset('WebPreferences')->search({user_id=>$user_id, page=>$page});
    $item->delete;
  }

sub initialize_basefile
  {

    my ($self, %opts) = self_or_default(@_);
    my $basename = $opts{basename};
    my $prog=$opts{prog} || "CoGe";
    my $return_name = $opts{return_name};
    my $tempdir = $opts{tempdir} || $TEMPDIR;
    if ($basename)
      {
#	print STDERR "Have basename: $basename\n";
	($basename) = $basename =~ /([^\/].*$)/;
	my ($x, $cleanname) = check_taint($basename);
	$self->basefilename($cleanname);
	$self->basefile($tempdir."/$prog/".$cleanname);
	$self->logfile($self->basefile.".log");
	$self->sqlitefile($self->basefile.".sqlite");
      }
    else
      {
	mkdir "$TEMPDIR/$prog",0777 unless -d "$TEMPDIR/$prog";
	my $file = new File::Temp ( TEMPLATE=>$prog.'_XXXXXXXX',
				    DIR=>"$TEMPDIR/$prog/",
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

1;
