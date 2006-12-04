# -*- perl -*-

# t/007_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Image' ); }

my $object = CoGeX::Image->new ();
isa_ok ($object, 'CoGeX::Image');


