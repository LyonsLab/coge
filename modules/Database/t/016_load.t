# -*- perl -*-

# t/016_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Feature' ); }

my $object = CoGeX::Feature->new ();
isa_ok ($object, 'CoGeX::Feature');


