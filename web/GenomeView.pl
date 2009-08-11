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
$ENV{PATH} = "/opt/apache2/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw( $PAGE_NAME $DATE $DEBUG $USER $FORM $coge);

$PAGE_NAME="GenomeView.pl";
$DEBUG = 0;
$FORM = new CGI;
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$coge = CoGeX->dbconnect();

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       grab_sequence=>\&grab_sequence,
		       save_options=>\&save_options,
		       get_genome_info=>\&get_genome_info,
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header; gen_html();



sub gen_html
  {
    my $html; #=  "Content-Type: text/html\n\n";
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(LOGO_PNG=>"GenomeView-logo.png");
    $template->param(TITLE=>'Genome Viewer');
    $template->param(PAGE_TITLE=>'Genome Viewer');
    $template->param(HELP=>'/wiki/index.php?title=GenomeView');
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $template->param(BOX_WIDTH=>"100%");
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
    my $chr = $form->param('chr') if $form->param('chr');
    my $dsid = $form->param('ds') if $form->param('ds');
    $dsid = $form->param('dsid') if $form->param('dsid');
    my $z = $form->param('z') if $form->param('z');
    my $loc = $form->param('x') if $form->param('x');
    my $ver = $form->param('ver') if $form->param('ver');
    my $org = $form->param('org') if $form->param('org');
    my $gstid = $form->param('gstid') if $form->param('gstid');
    my $fid = $form->param('fid') if $form->param('fid');
    my $dsgid = $form->param('dsgid') if $form->param('dsgid');

    my $prefs = load_settings(user=>$USER, page=>$PAGE_NAME);
    my $chr_length;

    if ($fid)
      {
	my $feat = $coge->resultset('Feature')->find($fid);
	$chr = $feat->chromosome;
	$dsid = $feat->dataset->id;
	$loc = $feat->start;
      }
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	$gstid = $dsg->type->id;
	my ($ds) = $dsg->datasets(chr=>$chr);
	$dsid = $ds->id;
      }
    if ($dsid)
      {
	my $dso = $coge->resultset('Dataset')->find($dsid);
	$org = $dso->organism->name;
	$ver = $dso->version;
	$chr_length = $dso->last_chromosome_position($chr);
	foreach my $dsg (sort {$a->genomic_sequence_type_id <=> $b->genomic_sequence_type_id} $dso->dataset_groups)
	  {
	    last if $dsgid;
	    $dsgid = $dsg->id if $gstid && $dsg->genomic_sequence_type_id == $gstid;
	    $dsgid = $dsg->id unless $gstid;
	  }
      }
    $loc = 1 unless $loc;
    $z=2 unless $z;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GenomeView.tmpl');
    $ver = "unk" unless $ver;
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
    my %default_true = (
			gc=>'true',
			genes=>'true',
		       );
    foreach my $item (qw (gc gaga gbox genes wobblegc wobble50gc localdup funcdomain prot))
      {
	my $show = $prefs->{$item};
	$show = $default_true{$item} unless $show;
	$show = 'false' unless $show;
	$template->param(uc($item)=>$show);
      }
    my $html = $template->output;
    my $org_name = "$org (v $ver), Chromosome: $chr, Dataset ID No. $dsid";
    return $html, $org_name;
  }

sub generate_box_name
{
  my $form = shift || $FORM;
  my $ds = $form->param('ds');
  my $chr = $form->param('chr');
  my $dso = $coge->resultset('Dataset')->find($ds);
  my $org = $dso->organism->name;
  my $ver = $dso->version;
  my $title = "$org (v $ver), Chromosome: $chr, Dataset ID No. $ds";
  return $title;
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
    my $item = save_settings(opts=>$opts, user=>$USER, page=>$PAGE_NAME);
  }

sub get_genome_info
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
#    print STDERR Dumper \%opts;
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
