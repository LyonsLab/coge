# -*- perl -*-

# t/012_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Sequence_type' ); }

my $object = CoGe::Genome::DB::Sequence_type->new ();
isa_ok ($object, 'CoGe::Genome::DB::Sequence_type');


