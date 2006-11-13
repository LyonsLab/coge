# -*- perl -*-

# t/005-5_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Dataset' ); }

my $object = CoGe::Genome::DB::Dataset->new ();

isa_ok ($object, 'CoGe::Genome::DB::Dataset');


