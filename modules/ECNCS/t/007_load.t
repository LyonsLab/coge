# -*- perl -*-

# t/007_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Classification' ); }

my $object = CoGe::ECNCS::DB::Classification->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Classification');
