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
$DIAGSDIR = "$DIR/diags";
$DOTPLOT = $DIR."/bin/dotplot.pl";

$FORM = new CGI;
$coge = CoGeX->dbconnect();

my $dsgid1 = $FORM->param('dsg1');
my $dsgid2 = $FORM->param('dsg2');
my $chr1 = $FORM->param('chr1');
my $chr2 = $FORM->param('chr2');
my $basename = $FORM->param('base');
my $url_only = $FORM->param('link');
my $square = $FORM->param('square');
my $flip = $FORM->param('flip');
my $regen = $FORM->param('regen');
my $width = $FORM->param('width') || 600;
my $grid = $FORM->param('grid');
my $ksdb = $FORM->param('ksdb');
my $kstype = $FORM->param('kstype');
my $log = $FORM->param('log');
my $min = $FORM->param('min');
my $max = $FORM->param('max');

$grid = 1 unless defined $grid;
$DEBUG=1 if $FORM->param('debug');
exit unless ($dsgid1 && $dsgid2 && $chr1 && $chr2 && $basename);
my ($md51, $md52, $mask1, $mask2, $type1, $type2, $blast,$params) = $basename =~/(.*?)_(.*?)\.(\d+)-(\d+)\.(\w+)-(\w+)\.(\w+)\.dag\.?a?l?l?_(.*)/ ;


my($dsg1) = $coge->resultset('DatasetGroup')->resolve($dsgid1);
my($dsg2) = $coge->resultset('DatasetGroup')->resolve($dsgid2);
exit unless $dsg1 && $dsg2;

my $name1 = $dsg1->organism->name;
$name1  =~ s/\///g;
$name1 =~ s/\s+/_/g;
$name1 =~ s/\(//g;
$name1 =~ s/\)//g;
$name1 =~ s/://g;
$name1 =~ s/;//g;
my $name2 = $dsg2->organism->name;
$name2  =~ s/\///g;
$name2 =~ s/\s+/_/g;
$name2 =~ s/\(//g;
$name2 =~ s/\)//g;
$name2 =~ s/://g;
$name2 =~ s/;//g;
my $dir = "$DIAGSDIR/$name1"."/".$name2;
my $dag_file = $dir."/".$basename;
#if ($ksdb)
#  {
#    ($ksdb) = $basename =~ /^(.*?CDS-CDS)/; 
#    $ksdb = $dir."/".$ksdb.".sqlite" if $ksdb; #that way it is not set if genomic sequence are compared.
#  }
$dag_file =~ s/\.dag_?.*//;
$dag_file .= ".dag.all";
my $outfile = $dir."/html/".$basename.".$chr1-$chr2.w$width";
my $res = generate_dotplot(dag=>$dag_file, coords=>"$dir/$basename.all.aligncoords", qchr=>$chr1, schr=>$chr2, 'outfile'=>$outfile, 'regen_images'=>'false', flip=>$flip, regen_images=>$regen, dsgid1=>$dsgid1, dsgid2=>$dsgid2, width=>$width, grid=>$grid, ksdb=>$ksdb, kstype=>$kstype, log=>$log, min=>$min, max=>$max);
if ($res)
  {
    $res=~s/$DIR/$URL/;
#    print STDERR $res,"\n";
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
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    my $qchr = $opts{qchr};
    my $schr = $opts{schr};
    my $flip = $opts{flip} || 0; 
    my $width = $opts{width} || 600;
    my $grid = $opts{grid} || 0; 
    $outfile .= ".flip" if $flip;
    my $regen_images = $opts{regen_images};
    my $ksdb = $opts{ksdb};
    my $kstype = $opts{kstype};
    my $log = $opts{log};
    my $min = $opts{min};
    my $max = $opts{max};
#    write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    my $cmd = $DOTPLOT;
    if ($ksdb && -r $ksdb)
      {
	$cmd .= qq{ -ksdb $ksdb -kst $kstype};
	$cmd .= qq{ -log 1} if $log;
	$cmd .= qq{ -min $min} if defined $min && $min =~/\d/;
	$cmd .= qq{ -max $max} if defined $max && $max =~/\d/;
	$outfile .= ".$kstype";
      }
    my $tmp = $outfile;
    $tmp .= ".$min" if defined $min && $min =~/\d/;
    $tmp .= ".$max" if defined $max && $max =~/\d/;
    if (-r $tmp.".html" && !$regen_images)
      {
#	write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }

    $cmd .= qq{ -d $dag -a $coords -b $outfile -l '' -dsg1 $dsgid1 -dsg2 $dsgid2 -w $width -lt 1 -chr1 $qchr -chr2 $schr -flip $flip -grid 1};
    print STDERR "Running: ",$cmd,"\n" if $DEBUG;
    `$cmd`;
    $outfile .= ".$min" if defined $min && $min =~/\d/;
    $outfile .= ".$max" if defined $max && $max =~/\d/;
    
    return $outfile if -r $outfile.".html";
  }
