# -*- perl -*-

# t/007_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Feature_name' ); }

my $object = CoGe::Genome::DB::Feature_name->new ();
isa_ok ($object, 'CoGe::Genome::DB::Feature_name');


