# -*- perl -*-

# t/010_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Location' ); }

my $object = CoGe::Genome::DB::Location->new ();
isa_ok ($object, 'CoGe::Genome::DB::Location');


