# -*- perl -*-

# t/004_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::AnnotationType' ); }

my $object = CoGeX::AnnotationType->new ();
isa_ok ($object, 'CoGeX::AnnotationType');


