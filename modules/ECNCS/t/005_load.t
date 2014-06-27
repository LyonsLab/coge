# -*- perl -*-

# t/005_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Algorithm_data' ); }

my $object = CoGe::ECNCS::DB::Algorithm_data->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Algorithm_data');
