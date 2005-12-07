# -*- perl -*-

# t/009_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Genomic_sequence' ); }

my $object = CoGe::Genome::DB::Genomic_sequence->new ();
isa_ok ($object, 'CoGe::Genome::DB::Genomic_sequence');


