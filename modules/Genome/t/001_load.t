# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB' ); }

my $object = CoGe::Genome::DB->new ();
isa_ok ($object, 'CoGe::Genome::DB');


