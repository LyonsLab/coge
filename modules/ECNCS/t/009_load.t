# -*- perl -*-

# t/009_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Status' ); }

my $object = CoGe::ECNCS::DB::Status->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Status');
