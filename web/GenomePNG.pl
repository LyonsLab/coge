#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Graphics;
use Data::Dumper;

$ENV{PATH} = "";
my $form = new CGI;
print STDERR $form->self_url(-full=>1),"\n";
my $gfx = new CoGe::Graphics;
my $start = $form->param('start') || $form->param('x') || $form->param('xmin') || 0;#28520458;
my $stop = $form->param('stop') || $form->param('xmax');# || 6948000;#6949600;#190000;
my $ds = $form->param('ds') || $form->param('di');
my $version = $form->param('v') || $form->param('version');
my $org = $form->param('o') || $form->param('org') ||$form->param('organism') || $form->param('org_id') || $form->param('oid');
my $chr = $form->param('chr') ||$form->param('chromosome');
my $iw = $form->param('iw') || $form->param('width') || $form->param('tile size')|| $form->param('tile_size') || 256;
my $ih = $form->param('ih') || $form->param('height');;
my $mag = $form->param('m') || $form->param('mag') || $form->param('magnification');
my $z = $form->param('z');
my $file = $form->param('file');# || "./tmp/pict.png";
my $start_pict = $form->param('start_pict');
my $simple = $form->param('simple');
my $chr_start_height = $form->param('csh') || 200;
my $chr_mag_height = $form->param('cmh') || 5;
my $feat_start_height = $form->param('fsh') || 10;
my $feat_mag_height = $form->param('fmh') || 2;
my $forcefit = $form->param('forcefit') || 0;
my $request = $form->param('request') || "get_image";
$forcefit = 1 if defined $start && defined $stop;
my @layers;
foreach my $layer ($form->param('layers'))
  {
    foreach my $item (split /,/, $layer)
      {
	push @layers, $item;
      }
  }
my @fids = $form->param('fid'); #used to highlight special features by their database id
my @fnames = $form->param('fn'); #used to highlight special features by their name
if ($request eq "get_annotation")
  {
    my $html;
    my $opts = "org=$org chr=$chr x=$start";
    $opts .= " version=$version" if $version;
    open (IN, "/opt/apache/CoGe/FeatAnno.pl $opts |");
    while (<IN>)
      {
        print $_;
      }
    close IN;
    
  }
else
  {
    CoGe::Graphics->genomic_view(
				 start=>$start,
				 stop=>$stop,
				 ds=>$ds,
				 version=>$version,
				 org=>$org,
				 chr=>$chr,
				 iw=>$iw,
				 ih=>$ih,
				 mag=>$mag,
				 z=>$z,
				 file=>$file,
				 start_pict=>$start_pict,
				 simple=>$simple,
				 csh=>$chr_start_height,
				 cmh=>$chr_mag_height,
				 fsh=>$feat_start_height,
				 fmh=>$feat_mag_height,
				 fids=>\@fids,
				 fns=>\@fnames,
				 forcefit=>$forcefit,
				 layers=>\@layers,
				);
  }
