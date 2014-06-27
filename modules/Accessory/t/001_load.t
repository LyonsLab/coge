# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Accessory::Tile::Cache' ); }

my $object = CoGe::Accessory::Tile::Cache->new ();
isa_ok ($object, 'CoGe::Accessory::Tile::Cache');
