# -*- perl -*-

# t/004_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Algorithm' ); }

my $object = CoGe::ECNCS::DB::Algorithm->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Algorithm');
