# -*- perl -*-

# t/022_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::FeatureListGroup' ); }

my $object = CoGeX::FeatureListGroup->new ();
isa_ok ($object, 'CoGeX::FeatureListGroup');


