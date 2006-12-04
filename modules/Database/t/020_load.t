# -*- perl -*-

# t/020_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::FeatureListConnector' ); }

my $object = CoGeX::FeatureListConnector->new ();
isa_ok ($object, 'CoGeX::FeatureListConnector');


