# -*- perl -*-

# t/014_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::Sequence' ); }

my $object = CoGeX::Sequence->new ();
isa_ok ($object, 'CoGeX::Sequence');


