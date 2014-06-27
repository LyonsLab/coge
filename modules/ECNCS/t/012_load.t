# -*- perl -*-

# t/012_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Ecncs' ); }

my $object = CoGe::ECNCS::DB::Ecncs->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Ecncs');
