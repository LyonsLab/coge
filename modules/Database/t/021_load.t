# -*- perl -*-

# t/021_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::UserGroupConnector' ); }

my $object = CoGeX::UserGroupConnector->new ();
isa_ok ($object, 'CoGeX::UserGroupConnector');


