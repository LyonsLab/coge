#! /usr/bin/perl -w

use strict;
use CoGe::Accessory::Tile::Cache;
use CGI;
use Data::Dumper;

my $form = new CGI;
my $img_generator_url = 'http://biocon.berkeley.edu/CoGe/GenomePNG.pl?';
my $x = $form->param('x') || $form->param('start');
my $z = $form->param('z');
my $iw = $form->param('iw') || $form->param('width') || $form->param('tile size')|| $form->param('tile_size');
my $di = $form->param('di');
my $chr = $form->param('chr');
my $org_id = $form->param('o') || $form->param('org') ||$form->param('organism') || $form->param('org_id') || $form->param('oid');
my $version = $form->param('v') || $form->param('version');
my $normalize = $form->param('normalize');
my $cache_file        = "/opt/apache/CoGe/cache/tile";

$cache_file .= ".iw";
$cache_file .= $iw if $iw;
$cache_file .= ".z";
$cache_file .= $z if $z;
$cache_file .= ".di";
$cache_file .= $di if $di;
$cache_file .= ".chr";
$cache_file .= $chr if $chr;
$cache_file .= ".oid";
$cache_file .= $org_id if $org_id;
$cache_file .= ".v";
$cache_file .= $version if $version;

my $base_units_per_pixel = 10/ $iw;

my $cache_object = CoGe::Accessory::Tile::Cache->new(
      $img_generator_url,$base_units_per_pixel,$cache_file
    );

# any parameters that should not affect the url used to generate the
# image are accessorized. all other parameters should go to the
# anonymous hash in get_tile()
$cache_object->DEBUG(1);

# to turn on the indexing by tile so &x=20
# finds the 20th tile
# see Tile::Cache for other accessors
$cache_object->normalize($normalize || 1);

print "Content-type: image/png\n\n";


# only parameters that are required to generate
# the image should go here

my %opts;
$opts{chr}   = $chr     if $chr;
$opts{di}    = $di      if $di;
$opts{org_id}= $org_id  if $org_id;
$opts{v}     = $version if $version;
print $cache_object->get_tile({
			       'x'         => $x,  # required
			       'z'         => $z,  # required
			       'tile_size' => $iw, # required
			       %opts,
			      });
