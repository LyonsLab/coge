# -*- perl -*-

# t/005-5_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Data_information' ); }

my $object = CoGe::Genome::DB::Data_information->new ();
isa_ok ($object, 'CoGe::Genome::DB::Data_information');


