#! /usr/bin/perl -w
use strict;
use CGI;
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
($USER) = CoGe::Accessory::LogUser->get_user();
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
print gen_html();

sub gen_html
  {
    my $html =  "Content-Type: text/html\n\n";
    unless ($USER)
      {
        $html .= login();
      }
    else
      {
        my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
#	$template->param(NO_JQUERY=>1);
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
    return $html;
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
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GenomeView.tmpl');
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

