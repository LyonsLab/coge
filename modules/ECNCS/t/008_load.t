# -*- perl -*-

# t/008_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Algorithm_run' ); }

my $object = CoGe::ECNCS::DB::Algorithm_run->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Algorithm_run');
