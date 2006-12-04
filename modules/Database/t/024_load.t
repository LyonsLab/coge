# -*- perl -*-

# t/024_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::FeatureListGroupImageConnector' ); }

my $object = CoGeX::FeatureListGroupImageConnector->new ();
isa_ok ($object, 'CoGeX::FeatureListGroupImageConnector');


