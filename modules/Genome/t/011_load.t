# -*- perl -*-

# t/011_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Sequence' ); }

my $object = CoGe::Genome::DB::Sequence->new ();
isa_ok ($object, 'CoGe::Genome::DB::Sequence');


