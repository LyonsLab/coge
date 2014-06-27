# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Algos::Pairwise' ); }

my $object = CoGe::Algos::Pairwise->new ();
isa_ok ($object, 'CoGe::Algos::Pairwise');
