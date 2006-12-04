# -*- perl -*-

# t/006_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::AnnotationTypeGroup' ); }

my $object = CoGeX::AnnotationTypeGroup->new ();
isa_ok ($object, 'CoGeX::AnnotationTypeGroup');


