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

use vars qw( $DATE $DEBUG $USER $FORM $coge);

$DEBUG = 0;
$FORM = new CGI;
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       grab_sequence=>\&grab_sequence,
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
        $template->param(USER=>$USER);
        $template->param(DATE=>$DATE);
        $template->param(BOX_NAME=>generate_box_name());
	$template->param(BODY_ONLOAD=>'init();');
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
    if ($ds)
      {
	my $dso = $coge->resultset('Dataset')->find($ds);
	$org = $dso->organism->name;
	$ver = $dso->version;
      }
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GeLo.tmpl');
    $template->param(CHR=>$chr);
    $template->param(VER=>$ver);
    $template->param(ORG=>$org);
    $template->param(DS=>$ds);
    $template->param(LOC=>$loc);
    $template->param(ZOOM=>$z);
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
  	print STDERR Dumper \@vals;
  	return $vals[0],$vals[1];
  }

