# -*- perl -*-

# t/006_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Feature' ); }

my $object = CoGe::Genome::DB::Feature->new ();
isa_ok ($object, 'CoGe::Genome::DB::Feature');


