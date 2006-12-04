# -*- perl -*-

# t/011_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Organism' ); }

my $object = CoGeX::Organism->new ();
isa_ok ($object, 'CoGeX::Organism');


