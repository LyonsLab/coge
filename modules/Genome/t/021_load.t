# -*- perl -*-

# t/015_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::Feature_list_connector' ); }

my $object = CoGe::Genome::DB::Feature_list_connector->new ();
isa_ok ($object, 'CoGe::Genome::DB::Feature_list_connector');


