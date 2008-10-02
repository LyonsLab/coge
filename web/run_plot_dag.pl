#!/usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use Data::Dumper;

$ENV{PATH} = "/opt/apache2/CoGe/";
umask(0);
use vars qw( $DATE $DEBUG $DIR $URL $USER $FORM $coge $cogeweb $DATADIR $DIAGSDIR $PYTHON $PLOT_DAG);

$DEBUG = 0;
$DIR = "/opt/apache/CoGe/";
$URL = "/CoGe/";
$DATADIR = "$DIR/data/";
$DIAGSDIR = "$DATADIR/diags";
$PYTHON = "/usr/bin/python";
$PLOT_DAG = $PYTHON ." ".$DIR."/bin/parepair/plot_dag.py";


$FORM = new CGI;
$coge = CoGeX->dbconnect();

my $org1 = $FORM->param('org1');
my $org2 = $FORM->param('org2');
my $chr1 = $FORM->param('chr1');
my $chr2 = $FORM->param('chr2');
my $basename = $FORM->param('base');
$DEBUG=1 if $FORM->param('debug');

exit unless ($org1 && $org2 && $chr1 && $chr2 && $basename);
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
	next unless $chr == $chr1;
	my $last = $ds->last_chromosome_position($chr);
	$params{1}= {dsid => $ds->id,
		     chr_end=>$last};
      }
  }

foreach my $ds (@ds2)
  {
    foreach my $chr ($ds->get_chromosomes)
      {
	next unless $chr == $chr2;
	my $last = $ds->last_chromosome_position($chr);
	$params{2}= {dsid => $ds->id,
		     chr_end=>$last};
      }
  }
my $name1 = $org1->name;
$name1 =~ s/\s+/_/g;
$name1 =~ s/\(//g;
$name1 =~ s/\)//g;

my $name2 = $org2->name;
$name2 =~ s/\s+/_/g;
$name2 =~ s/\(//g;
$name2 =~ s/\)//g;
my $dir = "$DIAGSDIR/$name1"."_".$name2;
my $dag_file = $dir."/".$basename;
$dag_file =~ s/\.dag_.*//;
$dag_file .= ".dag.all";
my $res = generate_dotplot(dag=>$dag_file, coords=>"$dir/$basename.all.aligncoords", qchr=>"a".$md51."_".$chr1, schr=>"b".$md52."_".$chr2, qdsid=>$params{1}{dsid}, q_max=>$params{1}{chr_end},  sdsid=>$params{2}{dsid}, s_max=>$params{2}{chr_end}, 'qlabel'=>$name1.":".$chr1, 'slabel'=>$name2.":".$chr2, 'outfile'=>$dir."/html/".$md51."_".$chr1."-".$mask1."-".$type1."_".$md52."_".$chr2."-".$mask2."-".$type2."_".$params.".".$blast.".html", 'regen_images'=>'false');

if ($res)
  {
    $res=~s/$DIR/$URL/;
    print qq{Content-Type: text/html


<html>
<head>
<meta HTTP-EQUIV="REFRESH" content="0; url=$res">
$res
</head>
</html>
}
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
    my $qchr = $opts{qchr};
    my $schr = $opts{schr};
    my $q_dsid = $opts{qdsid};
    my $s_dsid = $opts{sdsid};
    my $q_label = $opts{qlabel};
    my $s_label = $opts{slabel};
    my $q_max = $opts{"q_max"};
    my $s_max = $opts{"s_max"};
    my $regen_images = $opts{regen_images}=~/true/i ? 1 : 0;
#    print STDERR Dumper \%opts;
    if (-r $outfile && !$regen_images)
      {
#	write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }
    my $cmd = "$PLOT_DAG -d $dag -a $coords --html $outfile --q_dsid $q_dsid --s_dsid $s_dsid --qchr $qchr --schr $schr --q_max $q_max --s_max $s_max";
#    my $cmd = "$DAG_PLOT -d $dag -a $coords -f $outfile";
    $cmd .= " --q_label=$q_label" if $q_label;
    $cmd .= " --s_label=$s_label" if $s_label;
#    write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    `$cmd`;
    print STDERR $cmd,"\n" if $DEBUG;
    return $outfile if -r $outfile;
  }
