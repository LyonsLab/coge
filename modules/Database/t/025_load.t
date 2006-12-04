# -*- perl -*-

# t/025_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::UserSession' ); }

my $object = CoGeX::UserSession->new ();
isa_ok ($object, 'CoGeX::UserSession');


