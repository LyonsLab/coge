# -*- perl -*-

# t/010_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Author' ); }

my $object = CoGe::ECNCS::DB::Author->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Author');
