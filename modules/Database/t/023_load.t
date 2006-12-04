# -*- perl -*-

# t/023_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::UserGroupFeatureListPermissionConnector' ); }

my $object = CoGeX::UserGroupFeatureListPermissionConnector->new ();
isa_ok ($object, 'CoGeX::UserGroupFeatureListPermissionConnector');


