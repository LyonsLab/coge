#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGeX;
use CoGeX::Feature;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;

use vars qw( $TEMPDIR $USER $DATE $BASEFILE $coge $FORM);

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = "/opt/apache/CoGe/tmp/";
$FORM = new CGI;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
			);
$pj->js_encode_function('escape');
#print $pj->build_html($FORM, \&gen_html);
print $FORM->header;
print gen_html();

sub gen_html
  {
    my $html;    
    unless ($USER)
      {
	$html = login();
      }
    else
     {
    my ($body) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Feature List Viewer');
    $template->param(HELP=>'BLAST');
   # print STDERR "user is: ",$USER,"\n";
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(BOX_NAME=>'CoGe: Blast');
    $template->param(BODY=>$body);
    $html .= $template->output;
   }
 }
 
sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatList.tmpl');
    my $form = $FORM;
    $BASEFILE = $form->param('basename');
    my $table = generate_table();
    $template->param(INFO=>$table);
    $template->output;
  }
  
sub generate_table
  {
    my @table;
    my @feat_list;
    my $file = "$TEMPDIR/$BASEFILE.featlist";
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
    {
     push @feat_list,$_;
    }
    close IN;
    
    foreach my $info (@feat_list)
    {
      my ($featid,$dsid) = $info =~ /(\d+)\s+(\d+)/;
      print STDERR "($featid,$dsid)\n";
	 my ($ds) = $coge->resultset("Dataset")->find($dsid);
   	 my ($feat) = $coge->resultset("Feature")->find($featid);
   	 
   	 my ($name) = sort $feat->names;
   	 my $hpp = $feat->annotation_pretty_print_html();
   	 
   	 push @table,{NAME=>$name,TYPE=>$feat->type->name,LOC=>$feat->start."-".$feat->stop,CHR=>$feat->chr,STRAND=>$feat->strand,ORG=>$ds->organism->name."(v ".$feat->version.")",HPP=>$hpp};
    }
   return \@table;
  }