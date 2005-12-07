# -*- perl -*-

# t/013_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Organism' ); }

my $object = CoGe::Genome::DB::Organism->new ();
isa_ok ($object, 'CoGe::Genome::DB::Organism');


