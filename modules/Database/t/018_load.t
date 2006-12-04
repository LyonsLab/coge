# -*- perl -*-

# t/018_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::FeatureList' ); }

my $object = CoGeX::FeatureList->new ();
isa_ok ($object, 'CoGeX::FeatureList');


