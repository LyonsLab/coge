# -*- perl -*-

# t/003_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB::Spike' ); }

my $object = CoGe::ECNCS::DB::Spike->new ();
isa_ok ($object, 'CoGe::ECNCS::DB::Spike');
