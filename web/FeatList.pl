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
		       feats_parse=>\&feats_parse,
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();

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
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;#: $opts{feature_list};
    foreach my $item ($form->param('fid'))
      {
	push @$feat_list, $item if $item =~ /^\d+$/;
      }
    my $table = generate_table(feature_list=>$feat_list);
    if ($table)
      {
	$template->param(INFO=>$table);
	return $template->output;
      }
    else
      {
	return "No feature ids were specified.";
      }
  }
  
sub generate_table
  {
    my %opts = @_;
    my $feat_list = $opts{feature_list};
    return unless @$feat_list;
    my @table;
    my $count = 1;
    foreach my $info (@$feat_list)
    {
      my $featid = $info;
      #print STDERR $featid,"\n";
      my ($feat) = $coge->resultset("Feature")->find($featid);
      unless ($feat)
	{
	  warn "feature id $featid failed to return a valid feature object\n";
	  next;
	}
      my ($name) = sort $feat->names;
      my $hpp = $feat->annotation_pretty_print_html();
      my $row_style = $count%2 ? "even" : "odd";
      push @table,{FEATID=>$info,NAME=>$name,TYPE=>$feat->type->name,LOC=>$feat->start."-".$feat->stop,CHR=>$feat->chr,STRAND=>$feat->strand,ORG=>$feat->organism->name."(v ".$feat->version.")",HPP=>$hpp, TABLE_ROW=>$row_style};
      $count++;
    }
   return \@table;
  }
  
  sub read_file
  {
    my $file = "$TEMPDIR/$BASEFILE.featlist";
    my @featlist;
    unless (-r $file)
      {
	warn "unable to read file $file for feature ids\n";
	return \@featlist;
      }
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
	chomp;
	push @featlist,$_;
      }
    close IN;
    return \@featlist;
  }
  
  sub feats_parse #Send to GEvo
  {
    my $accn_list = shift;
    print STDERR $accn_list,"\n";
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my @list;
    my $url = "/CoGe/GEvo.pl?";
    my $count = 1;
    #print STDERR $url,"\n";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid$count=$featid&";
		$count ++;
      }
    $count--;
    return ("alert",$count) if $count > 10;
    $url .= "num_seqs=$count";
    return $url;
  }
