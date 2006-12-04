# -*- perl -*-

# t/019_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::UserGroup' ); }

my $object = CoGeX::UserGroup->new ();
isa_ok ($object, 'CoGeX::UserGroup');


