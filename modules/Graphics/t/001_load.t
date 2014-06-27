# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Graphics' ); }

my $object = CoGe::Graphics->new ();
isa_ok ($object, 'CoGe::Graphics');
