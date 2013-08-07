#!/usr/bin/perl -w

use strict;
use CGI;
use Data::Dumper;

use CoGe::Graphics;
use CoGe::Accessory::Web;
use CoGeX;
no warnings 'redefine';

$ENV{PATH} = "";
use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr);

$P = CoGe::Accessory::Web::get_defaults();

my $form = new CGI;

#print STDERR $form->self_url(-full=>1),"\n";
my $gfx = new CoGe::Graphics;
my $start =
     $form->param('start')
  || $form->param('x')
  || $form->param('xmin')
  || 0;    #28520458;
my $stop = $form->param('stop')
  || $form->param('xmax');    # || 6948000;#6949600;#190000;
my $ds      = $form->param('ds') || $form->param('di');
my $version = $form->param('v')  || $form->param('version');
my $org =
     $form->param('o')
  || $form->param('org')
  || $form->param('organism')
  || $form->param('org_id')
  || $form->param('oid');
my $chr = $form->param('chr') || $form->param('chromosome');
my $iw =
     $form->param('iw')
  || $form->param('width')
  || $form->param('tile size')
  || $form->param('tile_size')
  || 256;
my $ih     = $form->param('ih') || $form->param('height');
my $file   = $form->param('file');                        # || "./tmp/pict.png";
my $simple = $form->param('simple');
my $chr_height  = $form->param('ch') || 285;
my $feat_height = $form->param('fh') || 50;
my $request     = $form->param('request') || "get_image";
my $gstid       = $form->param('gstid');
my $dsgid       = $form->param('dsg');
my @ftid        = $form->param('ftid');
my @expid       = $form->param('expid');
my $red         = $form->param('r');
my $green       = $form->param('g');
my $blue        = $form->param('b');
my $color;
$color = [ $red, $green, $blue ]
  if defined $red && defined $green && defined $blue;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
my $coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

my @layers;
foreach my $layer ( $form->param('layers') ) {
    foreach my $item ( split /,/, $layer ) {
        push @layers, $item;
    }
}
$simple = 1 if @layers == 2 && $layers[0] =~ /background/i;
my @fids =
  $form->param('fid');  #used to highlight special features by their database id
my @fnames =
  $form->param('fn');    #used to highlight special features by their name

if ( $request eq "get_annotation" ) {
    my $html;
    my $opts = "org=$org chr=$chr x=$start";
    $opts .= " version=$version" if $version;
    open( IN, $P->{COGEDIR} . "FeatAnno.pl $opts |" );
    while (<IN>) {
        print $_;
    }
    close IN;

}
else {
    CoGe::Graphics->genomic_view(
        start             => $start,
        stop              => $stop - 1,
        ds                => $ds,
        gstid             => $gstid,
        dsgid             => $dsgid,
        ftid              => \@ftid,
        expid             => \@expid,
        color             => $color,
        version           => $version,
        org               => $org,
        chr               => $chr,
        iw                => $iw,
        ih                => $ih,
        file              => $file,
        simple            => $simple,
        ch                => $chr_height,
        fh                => $feat_height,
        fids              => \@fids,
        fns               => \@fnames,
        layers            => \@layers,
        DEBUG             => 0,
        major_tick_labels => 1,
        minor_tick_labels => -1,
        coge              => $coge,
    );
}
