# -*- perl -*-

# t/013_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Permission' ); }

my $object = CoGeX::Permission->new ();
isa_ok ($object, 'CoGeX::Permission');


