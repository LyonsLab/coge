# -*- perl -*-

# t/002_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Annotation' ); }

my $object = CoGeX::Annotation->new ();
isa_ok ($object, 'CoGeX::Annotation');


