# -*- perl -*-

# t/002_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::ECNCS::DB' ); }

my $object = CoGe::ECNCS::DB->new ();
isa_ok ($object, 'CoGe::ECNCS::DB');
