# -*- perl -*-

# t/009_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Location' ); }

my $object = CoGeX::Location->new ();
isa_ok ($object, 'CoGeX::Location');


