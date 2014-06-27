# -*- perl -*-

# t/011_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Location' ); }

my $object = CoGe::ECNCS::DB::Location->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Location');
