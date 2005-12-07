# -*- perl -*-

# t/005_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Data_source' ); }

my $object = CoGe::Genome::DB::Data_source->new ();
isa_ok ($object, 'CoGe::Genome::DB::Data_source');


