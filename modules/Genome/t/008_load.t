# -*- perl -*-

# t/008_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Feature_type' ); }

my $object = CoGe::Genome::DB::Feature_type->new ();
isa_ok ($object, 'CoGe::Genome::DB::Feature_type');


