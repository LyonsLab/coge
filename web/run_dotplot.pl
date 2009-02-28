#!/usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use Data::Dumper;

$ENV{PATH} = "/opt/apache2/CoGe/";
umask(0);
use vars qw( $DATE $DEBUG $DIR $URL $USER $FORM $coge $cogeweb $DATADIR $DIAGSDIR $DOTPLOT);

$DEBUG = 0;
$DIR = "/opt/apache/CoGe/";
$URL = "/CoGe/";
$DATADIR = "$DIR/data/";
$DIAGSDIR = "$DATADIR/diags";
$DOTPLOT = $DIR."/bin/dotplot.pl";

$FORM = new CGI;
$coge = CoGeX->dbconnect();

my $org1 = $FORM->param('org1');
my $org2 = $FORM->param('org2');
my $chr1 = $FORM->param('chr1');
my $chr2 = $FORM->param('chr2');
my $basename = $FORM->param('base');
my $url_only = $FORM->param('link');
my $square = $FORM->param('square');
my $flip = $FORM->param('flip');
my $regen = $FORM->param('regen');
my $width = $FORM->param('width') || 600;
my $grid = $FORM->param('grid');
$grid = 1 unless defined $grid;
$DEBUG=1 if $FORM->param('debug');
unless ($org1 && $org2 && $chr1 && $chr2 && $basename)
  {
    warn "Problem running run_dotplot.pl\n";
  }
my ($md51, $md52, $mask1, $mask2, $type1, $type2, $blast,$params) = $basename =~/(.*?)_(.*?)\.(\d+)-(\d+)\.(\w+)-(\w+)\.(\w+)\.dag_(.*)/ ;


($org1) = $coge->resultset('Organism')->resolve($org1);
my @ds1 = $org1->current_datasets(type=>$mask1);
($org2) = $coge->resultset('Organism')->resolve($org2);
my @ds2 = $org2->current_datasets(type=>$mask2);

my %params;



foreach my $ds (@ds1)
  {
    foreach my $chr ($ds->get_chromosomes)
      {
	my @tmp = split/_/,$chr; #this stuff is a hack to deal with chromosome names liek supercontig_1 which are passed in as "1" . . .
	next unless $chr eq $chr1 || $chr1 eq $tmp[-1];
	$chr1 = $chr if $chr ne $chr1;
	my $last = $ds->last_chromosome_position($chr);
	$params{1}= {dsid => $ds->id,
		     chr_end=>$last};
      }
  }

foreach my $ds (@ds2)
  {
    foreach my $chr ($ds->get_chromosomes)
      {
	my @tmp = split/_/,$chr;
	next unless $chr eq $chr2 || $chr2 eq $tmp[-1];
	$chr2 = $chr if $chr ne $chr2;
	my $last = $ds->last_chromosome_position($chr);
	$params{2}= {dsid => $ds->id,
		     chr_end=>$last};
      }
  }

my $name1 = $org1->name;
$name1 =~ s/\s+/_/g;
$name1 =~ s/\(//g;
$name1 =~ s/\)//g;
$name1 =~ s/://g;
$name1 =~ s/\///g;
my $name2 = $org2->name;
$name2 =~ s/\s+/_/g;
$name2 =~ s/\(//g;
$name2 =~ s/\)//g;
$name2 =~ s/://g;
$name2 =~ s/\///g;
my $dir = "$DIAGSDIR/$name1"."/".$name2;
my $dag_file = $dir."/".$basename;
$dag_file =~ s/\.dag_.*//;
$dag_file .= ".dag.all";
my $outfile = $dir."/html/".$md51."_".$chr1."-".$mask1."-".$type1."_".$md52."_".$chr2."-".$mask2."-".$type2."_".$params.".".$blast.".w$width";
my $res = generate_dotplot(dag=>$dag_file, coords=>"$dir/$basename.all.aligncoords", qchr=>$chr1, schr=>$chr2, 'outfile'=>$outfile, 'regen_images'=>'false', flip=>$flip, regen_images=>$regen, oid1=>$org1->id, oid2=>$org2->id, width=>$width, grid=>$grid);
if ($res)
  {
    $res=~s/$DIR/$URL/;
    print $res if $url_only;
    print qq{Content-Type: text/html


<html>
<head>
<meta HTTP-EQUIV="REFRESH" content="0; url=$res.html">
$res
</head>
</html>
};
    #print $file;
  }

#input file
#5883268adf8ee5509a1b20e8bac20bb3_5883268adf8ee5509a1b20e8bac20bb3.2-2.CDS-CDS.blastn.dag_D400000_g200000_A3.all.aligncoords

#outputfile
#5883268adf8ee5509a1b20e8bac20bb3_9_5883268adf8ee5509a1b20e8bac20bb3_9_D400000_g200000_A3.blastn.html




sub generate_dotplot
  {
    my %opts = @_;
    my $dag = $opts{dag};
    my $coords = $opts{coords};
    my $outfile = $opts{outfile};
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    my $qchr = $opts{qchr};
    my $schr = $opts{schr};
    my $flip = $opts{flip} || 0; 
    my $width = $opts{width} || 600;
    my $grid = $opts{grid} || 0; 
    $outfile .= ".flip" if $flip;
    my $regen_images = $opts{regen_images};
    if (-r $outfile.".html" && !$regen_images)
      {
#	write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }
#    write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    my $cmd = qq{$DOTPLOT -d $dag -a $coords -b $outfile -l '' -o1 $oid1 -o2 $oid2 -w $width -lt 1 -chr1 $qchr -chr2 $schr -flip $flip -grid 1};
    `$cmd`;
    return $outfile if -r $outfile.".html";
  }
