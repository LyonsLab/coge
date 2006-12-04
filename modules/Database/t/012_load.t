# -*- perl -*-

# t/012_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Dataset' ); }

my $object = CoGeX::Dataset->new ();
isa_ok ($object, 'CoGeX::Dataset');


