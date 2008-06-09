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
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print gen_html();



sub gen_html
  {
    my $html; #=  "Content-Type: text/html\n\n";
    unless ($USER)
      {
        $html .= login();
      }
    else
      {
        my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
#	$template->param(HEADER_LINK=>'/CoGe/GeLo.pl');
        $template->param(LOGO_PNG=>"GeLo-logo.png");
        $template->param(TITLE=>'Genome Viewer');
        $template->param(HELP=>'');
        my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

	$template->param(LOGON=>1) unless $USER->user_name eq "public";
        $template->param(DATE=>$DATE);
        $template->param(BOX_NAME=>generate_box_name());
#	$template->param(BODY_ONLOAD=>'init();');
        $template->param(BODY=>gen_body());
        $html .= $template->output;
      }
    #return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $chr = $form->param('chr');
    my $ds = $form->param('ds');
    my $z = $form->param('z');
    my $loc = $form->param('x');
    my $ver = $form->param('ver');
    my $org = $form->param('org');
    my $prefs = load_settings(user=>$USER, page=>$PAGE_NAME);
    if ($ds)
      {
	my $dso = $coge->resultset('Dataset')->find($ds);
	$org = $dso->organism->name;
	$ver = $dso->version;
      }
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GeLo.tmpl');
    $ver = "unk" unless $ver;
    $template->param(CHR=>$chr);
    $template->param(VER=>$ver);
    $template->param(ORG=>$org);
    $template->param(DS=>$ds);
    $template->param(LOC=>$loc);
    $template->param(ZOOM=>$z);
    $template->param(SAVE_SETTINGS=>1) unless $USER->user_name eq "public";
    $template->param(FLAT=>"checked") if $prefs->{'flat'} eq "true";
    $template->param(EXPAND=>"checked") if $prefs->{'expand'} eq "true";
    $template->param(POPUPANNO=>"checked") if $prefs->{'popupanno'} eq "true";
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
    return $html;
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
