#!/usr/bin/perl -w

use strict;
use CGI;
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::Web;
no warnings 'redefine';

umask(0);
use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $DIR $URL $USER $FORM $coge $cogeweb $DATADIR $DIAGSDIR $DOTPLOT);
$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$DEBUG     = 0;
$DIR       = $P->{COGEDIR};
$URL       = $P->{URL};
$DATADIR   = $P->{DATADIR};
$DIAGSDIR  = $P->{DIAGSDIR};
$DOTPLOT   = $P->{DOTPLOT} . " -cf " . $ENV{COGE_HOME} . 'coge.conf';

$FORM   = new CGI;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

my $dsgid1       = $FORM->param('dsg1');
my $dsgid2       = $FORM->param('dsg2');
my $chr1         = $FORM->param('chr1');
my $chr2         = $FORM->param('chr2');
my $basename     = $FORM->param('base');
my $url_only     = $FORM->param('link');
my $square       = $FORM->param('square');
my $flip         = $FORM->param('flip');
my $regen        = $FORM->param('regen');
my $width        = $FORM->param('width') || 600;
my $grid         = $FORM->param('grid');
my $ksdb         = $FORM->param('ksdb');
my $kstype       = $FORM->param('kstype');
my $log          = $FORM->param('log');
my $min          = $FORM->param('min');
my $max          = $FORM->param('max');
my $metric       = $FORM->param('am');
my $relationship = $FORM->param('ar');
my $box_diags    = $FORM->param('bd');
my $color_type   = $FORM->param('ct');
my $color_scheme = $FORM->param('cs');
my $fid1         = $FORM->param('fid1');
my $fid2         = $FORM->param('fid2');
$box_diags = 1 if $box_diags && $box_diags eq "true";
$grid = 1 unless defined $grid;
$DEBUG = 1 if $FORM->param('debug');
exit
  unless ( $dsgid1 && $dsgid2 && defined $chr1 && defined $chr2 && $basename );
my ( $md51, $md52, $mask1, $mask2, $type1, $type2, $blast, $params ) =
  $basename =~
  /(.*?)_(.*?)\.(\d+)-(\d+)\.(\w+)-(\w+)\.(\w+)\.dag\.?a?l?l?_(.*)/;

my ($dsg1) = $coge->resultset('Genome')->resolve($dsgid1);
my ($dsg2) = $coge->resultset('Genome')->resolve($dsgid2);
exit unless $dsg1 && $dsg2;

my $name1 = $dsg1->organism->name;
$name1 =~ s/\///g;
$name1 =~ s/\s+/_/g;
$name1 =~ s/\(//g;
$name1 =~ s/\)//g;
$name1 =~ s/://g;
$name1 =~ s/;//g;
$name1 =~ s/'//g;
$name1 =~ s/"//g;
my $name2 = $dsg2->organism->name;
$name2 =~ s/\///g;
$name2 =~ s/\s+/_/g;
$name2 =~ s/\(//g;
$name2 =~ s/\)//g;
$name2 =~ s/://g;
$name2 =~ s/;//g;
$name2 =~ s/'//g;
$name2 =~ s/"//g;
my $dir      = "$DIAGSDIR/$name1" . "/" . $name2;
my $dag_file = $dir . "/" . $basename;

#if ($ksdb)
#  {
#    ($ksdb) = $basename =~ /^(.*?CDS-CDS)/;
#    $ksdb = $dir."/".$ksdb.".sqlite" if $ksdb; #that way it is not set if genomic sequence are compared.
#  }
$dag_file =~ s/\.dag_?.*//;
$dag_file .= ".dag.all";
my $outfile;
$outfile = ".$chr1-$chr2.w$width";
$outfile =~ s/\///g;
$outfile =~ s/\s+/_/g;
$outfile =~ s/\(//g;
$outfile =~ s/\)//g;
$outfile =~ s/://g;
$outfile =~ s/;//g;
$outfile = $dir . "/html/" . $basename . $outfile;

#my $res = generate_dotplot(dag=>$dag_file, coords=>"$dir/$basename.all.aligncoords", qchr=>$chr1, schr=>$chr2, 'outfile'=>$outfile, 'regen_images'=>'false', flip=>$flip, regen_images=>$regen, dsgid1=>$dsgid1, dsgid2=>$dsgid2, width=>$width, grid=>$grid, ksdb=>$ksdb, kstype=>$kstype, log=>$log, min=>$min, max=>$max, metric=>$metric);
my $res = generate_dotplot(
    dag            => $dag_file,
    coords         => "$dir/$basename",
    qchr           => $chr1,
    schr           => $chr2,
    'outfile'      => $outfile,
    'regen_images' => 'false',
    flip           => $flip,
    regen_images   => $regen,
    dsgid1         => $dsgid1,
    dsgid2         => $dsgid2,
    width          => $width,
    grid           => $grid,
    ksdb           => $ksdb,
    kstype         => $kstype,
    log            => $log,
    min            => $min,
    max            => $max,
    metric         => $metric,
    relationship   => $relationship,
    box_diags      => $box_diags,
    color_type     => $color_type,
    color_scheme   => $color_scheme,
    fid1           => $fid1,
    fid2           => $fid2
);
if ($res) {
    $res =~ s/$DIR/$URL/;

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

sub generate_dotplot {
    my %opts    = @_;
    my $dag     = $opts{dag};
    my $coords  = $opts{coords};
    my $outfile = $opts{outfile};
    my $dsgid1  = $opts{dsgid1};
    my $dsgid2  = $opts{dsgid2};
    my $fid1    = $opts{fid1};
    my $fid2    = $opts{fid2};
    my $qchr    = $opts{qchr};
    my $schr    = $opts{schr};
    my $flip    = $opts{flip} || 0;
    my $width   = $opts{width} || 600;
    my $grid    = $opts{grid} || 0;
    $outfile .= ".flip" if $flip;
    my $regen_images = $opts{regen_images};
    my $ksdb         = $opts{ksdb};
    my $kstype       = $opts{kstype};
    my $log          = $opts{log};
    my $min          = $opts{min};
    my $max          = $opts{max};
    my $metric       = $opts{metric};
    my $relationship = $opts{relationship};
    my $box_diags    = $opts{box_diags};
    my $color_type   = $opts{color_type};
    my $color_scheme = $opts{color_scheme};

#   CoGe::Accessory::Web::write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    my $cmd = $DOTPLOT;
    if ( $ksdb && -r $ksdb ) {
        $cmd     .= qq{ -ksdb $ksdb -kst $kstype};
        $cmd     .= qq{ -log 1} if $log;
        $cmd     .= qq{ -min $min} if defined $min && $min =~ /\d/;
        $cmd     .= qq{ -max $max} if defined $max && $max =~ /\d/;
        $outfile .= ".$kstype";
    }
    $outfile .= ".gene" if $metric       =~ /gene/i;
    $outfile .= ".s"    if $relationship =~ /s/i;
    $outfile .= ".box"  if $box_diags;
    $outfile .= ".ct$color_type"   if $color_type;
    $outfile .= ".cs$color_scheme" if defined $color_scheme;
    $outfile .= ".$min"            if defined $min && $min =~ /\d/;
    $outfile .= ".$max"            if defined $max && $max =~ /\d/;
    $outfile .= ".$fid1"           if $fid1;
    $outfile .= ".$fid2"           if $fid2;

    my $tmp = $outfile;
    $tmp .= ".$min" if defined $min && $min =~ /\d/;
    $tmp .= ".$max" if defined $max && $max =~ /\d/;
    if ( -r $tmp . ".html" && !$regen_images ) {

#CoGe::Accessory::Web::write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
        return $outfile;
    }

    $cmd .=
qq{ -d $dag -a "$coords" -b "$outfile" -dsg1 $dsgid1 -dsg2 $dsgid2 -w $width -lt 1 -chr1 "$qchr" -chr2 "$schr" -flip $flip -grid 1};
    $cmd .= qq{ -fid1 $fid1}       if $fid1;
    $cmd .= qq{ -fid2 $fid2}       if $fid2;
    $cmd .= qq { -am $metric}      if $metric;
    $cmd .= qq { -fb}              if $relationship && $relationship =~ /s/i;
    $cmd .= qq { -cdt $color_type} if $color_type;
    $cmd .= qq{ -bd 1}             if $box_diags;
    $cmd .= qq{ -color_scheme $color_scheme} if defined $color_scheme;
    print STDERR "Running: ", $cmd, "\n" if $DEBUG;
    my $x;
    ( $x, $cmd ) = CoGe::Accessory::Web::check_taint($cmd);

    #($cmd) = $cmd =~ /(.*)/;
    `$cmd` if $cmd;
    return $outfile if -r $outfile . ".html";
}
