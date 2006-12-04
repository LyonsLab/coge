# -*- perl -*-

# t/003_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::FeatureType' ); }

my $object = CoGeX::FeatureType->new ();
isa_ok ($object, 'CoGeX::FeatureType');


