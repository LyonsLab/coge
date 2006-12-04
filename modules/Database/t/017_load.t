# -*- perl -*-

# t/017_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::User' ); }

my $object = CoGeX::User->new ();
isa_ok ($object, 'CoGeX::User');


