# -*- perl -*-

# t/013_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::Accessory::Annotation' ); }

my $object = CoGe::Genome::Accessory::Annotation->new ();
isa_ok ($object, 'CoGe::Genome::Accessory::Annotation');


