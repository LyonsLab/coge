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
use Spreadsheet::WriteExcel;

use vars qw( $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM);

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
		       export_fasta_file=>\&export_fasta_file,
		       generate_excel_file=>\&generate_excel_file,
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
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;#: $opts{feature_list};
    foreach my $item ($form->param('fid'))
      {
	push @$feat_list, $item if $item =~ /^\d+$/;
      }
    my $dsid = $form->param('dsid') if $form->param('dsid');
    my $ftid = $form->param('ftid') if $form->param('ftid');
    push @$feat_list, @{get_fids_from_dataset(dsid=>$dsid, ftid=>$ftid)} if $dsid;
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

sub get_fids_from_dataset
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $ftid = $opts{ftid};
    my $search = {dataset_id=>$dsid};
    my $join={};
    if ($ftid)
      {
	$search->{feature_type_id}=$ftid;
      }
    my @ids = map{$_->id}$coge->resultset('Feature')->search($search);
    return \@ids;
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
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
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
  
  sub export_fasta_file
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "FastaView.pl?";
    foreach my $featid (split /,/,$accn_list)
    {
		$url .= "featid=$featid&";
	}
	$url =~s/&$//;
	return $url;
  }
  
sub generate_excel_file
  {
  	my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
  	$cogeweb = initialize_basefile(prog=>"FeatList");
  	my $basename = $cogeweb->basefile;
  	my ($filename) = $basename =~ /^\/opt\/apache\/CoGe\/tmp\/(FeatList_\w+)$/;
  	print STDERR $filename,"\n";
  	my $workbook = Spreadsheet::WriteExcel->new("$TEMPDIR/Excel_$filename.xls");
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i = 1;
       	 
   	 $worksheet->write(0,0,"Feature Name");
   	 $worksheet->write(0,1,"Type");
   	 $worksheet->write(0,2,"Location");
   	 $worksheet->write(0,3,"Strand");
   	 $worksheet->write(0,4,"Chromosome");
   	 $worksheet->write(0,5,"Organism (version)");
   	 $worksheet->write(0,6,"More information");
   	
   	foreach my $featid (split /,/,$accn_list)
    {
   	   my ($feat) = $coge->resultset("Feature")->find($featid);
   	   next unless $feat;
   	   my ($name) = sort $feat->names;
   	   my $app = $feat->annotation_pretty_print();
   	   $worksheet->write($i,0,"http://toxic.berkeley.edu/CoGe/FeatView.pl?accn=$name",$name);
   	   $worksheet->write($i,1,$feat->type->name);
   	   $worksheet->write($i,2,$feat->start."-".$feat->stop);
   	   $worksheet->write($i,3,$feat->strand);
   	   $worksheet->write($i,4,$feat->chr);
   	   $worksheet->write($i,5,$feat->organism->name."(v ".$feat->version.")");
   	   $worksheet->write($i,6,$app);
   	   $i++;
   	}
   	
   	$workbook->close() or die "Error closing file: $!";
   	#print STDERR "tmp/Excel_$filename.xls\n";
   	return "tmp/Excel_$filename.xls";
  }
