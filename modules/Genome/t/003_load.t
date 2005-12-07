# -*- perl -*-

# t/003_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Annotation_type' ); }

my $object = CoGe::Genome::DB::Annotation_type->new ();
isa_ok ($object, 'CoGe::Genome::DB::Annotation_type');


