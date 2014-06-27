# -*- perl -*-

# t/006_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Data_mask' ); }

my $object = CoGe::ECNCS::DB::Data_mask->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Data_mask');
