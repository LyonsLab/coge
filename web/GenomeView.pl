#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CGI::Carp 'fatalsToBrowser';
use HTML::Template;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use DBI;
use Data::Dumper;
no warnings 'redefine';

delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $DATE $DEBUG $USER $FORM $coge $COOKIE_NAME);
$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');
$ENV{PATH} = $P->{COGEDIR};
$PAGE_NAME="GenomeView.pl";
$DEBUG = 0;
$FORM = new CGI;
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));


$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) =$FORM->param('ticket');
CoGe::Accessory::Web->login_cas(ticket=>$cas_ticket, coge=>$coge, this_url=>$FORM->url()) if($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name=>$COOKIE_NAME,coge=>$coge);
my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       grab_sequence=>\&grab_sequence,
		       save_options=>\&save_options,
		       get_genome_info=>\&get_genome_info,
		      );
#$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header; gen_html();

sub gen_html
  {
    my $html; #=  "Content-Type: text/html\n\n";
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
    $template->param(LOGO_PNG=>"GenomeView-logo.png");
#    $template->param(TITLE=>'Genome Viewer');
    $template->param(PAGE_TITLE=>'Genome Viewer');
    $template->param(HELP=>'/wiki/index.php?title=GenomeView');
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
#    $template->param(BOX_WIDTH=>"100%");
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    my ($body, $org_name) = gen_body();
    $template->param(BODY=>$body);
    $template->param(BOX_NAME=>$org_name);
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my ($chr, $dsid, $z, $loc, $gstid, $fid, $dsgid, $show_legend, $prefs);
    $chr = $form->param('chr') if $form->param('chr');
    $dsid = $form->param('ds') if $form->param('ds');
    $dsid = $form->param('dsid') if $form->param('dsid');
    $z = $form->param('z') if defined ($form->param('z'));
    $loc = $form->param('x') if $form->param('x');
    $loc =~ s/,|\.//g if $loc; #remove commas and points if present
    $gstid = $form->param('gstid') if $form->param('gstid');
    $fid = $form->param('fid') if $form->param('fid');
    $dsgid = $form->param('dsgid') if $form->param('dsgid');
    $show_legend = $form->param('sl') if $form->param('sl');
    $prefs =CoGe::Accessory::Web::load_settings(user=>$USER, page=>$PAGE_NAME);
    my ($ds, $dsg, $gst);
    if ($fid)
      {
	my $feat = $coge->resultset('Feature')->find($fid);
	$chr = $feat->chromosome;
	$ds = $feat->dataset;
	$loc = $feat->start;
      }
    if ($dsgid)
      {
	$dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	$gst = $dsg->type;
	($ds) = $dsg->datasets(chr=>$chr);
	my @ds = $dsg->datasets();
      }
    
    if ($dsid)
      {
	$ds = $coge->resultset('Dataset')->find($dsid);
      }
    unless ($dsg)
      {
	foreach my $dsgt (sort {$a->genomic_sequence_type_id <=> $b->genomic_sequence_type_id} $ds->dataset_groups)
	  {
	    last if $dsgid;
	    if ($gstid && $dsgt->genomic_sequence_type_id == $gstid)
	      {
		$dsg = $dsgt;
		last;
	      }
	    $dsg = $dsgt;
	  }
      }
    $gst = $dsg->type;
    $dsid = $ds->id;
    $dsgid = $dsg->id;
    $gstid = $gst->id;
    my $ver = $dsg->version;
    my $org = $ds->organism->name;
    my $chr_length = $ds->last_chromosome_position($chr);

    my @feat_types;
    my $query = qq{select distinct(feature_type_id) from feature where dataset_id = $dsid};
    my $dbh = DBI->connect($connstr,$DBUSER,$DBPASS);
    my $sth = $dbh->prepare($query);
    $sth->execute;
    while (my $row = $sth->fetchrow_arrayref)
	  {
	    my $ftid = $row->[0];
	    push @feat_types, $coge->resultset('FeatureType')->find($ftid);
	  }
    $loc = 1 unless $loc;
    $z=5 unless defined $z;
    $z = 0 if $z < 0;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'GenomeView.tmpl');
    #set tiler program
    $template->param(TILE_SERVER=>$P->{SERVER}."tiler.pl?");
    #set layers
    foreach my $ft (@feat_types)
      {
	if ($ft->name eq "gene_space")
	  {
	    $template->param(GENE_SPACE_LAYER=>1);
	  }
	elsif ($ft->name =~ /duplicate/i)
	  {
	    $template->param(LOCAL_DUP_LAYER=>1);
	  }
	elsif ($ft->name =~ /domain/i)
	  {
	    $template->param(FUNC_DOMAIN_LAYER=>1);
	  }
	elsif ($ft->name =~ /repeat/i)
	  {
	    $template->param(REPEATS_LAYER=>1);
	  }
	elsif ($ft->name =~ /transposable/i)
	  {
	    $template->param(TE_LAYER=>1);
	  }
      }
    my ($gevo_group) = $coge->resultset('AnnotationTypeGroup')->search({name=>"gevo link"});
    if ($gevo_group)
      {
	my ($anno) = $coge->resultset('Annotation')->count({'feature.dataset_id'=>$dsid,'annotation_type.annotation_type_group_id'=>$gevo_group->id},{join=>['feature','annotation_type'],limit=>1});
	$template->param(GEVO_LINK_LAYER=>1) if ($anno);
      }


    my ($tandem_group) = $coge->resultset('AnnotationTypeGroup')->search({name=>"Tandem duplicates"});
    if ($tandem_group)
      {
	my ($anno) = $coge->resultset('Annotation')->count({'feature.dataset_id'=>$dsid,'annotation_type.annotation_type_group_id'=>$tandem_group->id},{join=>['feature','annotation_type'], limit=>1});
	$template->param(LOCAL_DUP_LAYER=>1) if $anno;
      }
    
    $template->param(CHR=>$chr);
    $template->param(VER=>$ver);
    $template->param(ORG=>$org);
    $template->param(DS=>$dsid);
    $template->param(DSG=>$dsgid);
    $template->param(LOC=>$loc);
    $template->param(ZOOM=>$z);
    $template->param(GSTID=>$gstid);
    $template->param(CHR_LENGTH=>$chr_length);
    $template->param(SAVE_SETTINGS=>1) unless $USER->user_name eq "public";
    $template->param(FLAT=>"checked") if $prefs->{'flat'}  && $prefs->{'flat'} eq "true";
    $template->param(EXPAND=>"checked") if $prefs->{'expand'}  && $prefs->{'expand'} eq "true";
    $template->param(POPUPANNO=>"checked") if $prefs->{'popupanno'}  && $prefs->{'popupanno'} eq "true";
    $template->param(SHOW_LEGEND=>1) if $show_legend;
    #stuff for people coming from maizegdb
    if ($ENV{HTTP_REFERER} && $ENV{HTTP_REFERER}=~/maizegdb/)
      {
	$template->param(MAIZEGDB=>1);
	$template->param(EXPAND=>"checked");
      }
    my %default_true = (
			gc=>'true',
			genes=>'true',
		       );
    foreach my $item (qw (gc gaga gbox genes wobblegc wobble50gc localdup funcdomain prot repeats other gevo_link))
      {
	my $show = $prefs->{$item};
	$show = $default_true{$item} unless $show;
	$show = 'false' unless $show;
	$template->param(uc($item)=>$show);
      }
    my $html = $template->output;
    my $org_name = "<span class=link onclick=window.open('OrganismView.pl?dsgid=$dsgid')>$org (v$ver),";
    $org_name .= " ".$dsg->name if $dsg->name;
    $org_name .= ": ".$dsg->description if $dsg->description;
    $org_name.= " Chromosome: $chr ".$gst->name." (dsgid$dsgid dsid$dsid)</span>";
    return $html, $org_name;
  }

sub grab_sequence
  {
  	my $new_value = shift;
  	my $first_value = shift;
  	my @vals = ($new_value,$first_value);
  	@vals = sort {$a <=> $b} (@vals);
  	return $vals[0],$vals[1];
  }

sub save_options
  {
    my %opts = @_;
    my $opts = Dumper \%opts;
    my $item =CoGe::Accessory::Web::save_settings(opts=>$opts, user=>$USER, page=>$PAGE_NAME);
  }

sub get_genome_info
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    return " " unless $dsgid;
    my $dsg = $coge->resultset("DatasetGroup")->find($dsgid);
    return "Unable to create dataset_group object for id: $dsgid" unless $dsg;
    my $html = "<span class=species>Datasets:</span><br><table>";
    $html .= "<tr><th>ID<th>Name<th>Description";
    foreach my $ds ($dsg->datasets)
      {
	$html .="<tr><td>".$ds->id."<td><a href=OrganismView.pl?dsid=".$ds->id." target=_new>".$ds->name."</a><td>".$ds->description;
      }
    $html .= "</table>";
    return $html;
  }
sub commify {
        my $input = shift;
        $input = reverse $input;
        $input =~ s<(\d\d\d)(?=\d)(?!\d*\.)><$1,>g;
        return scalar reverse $input;
}
