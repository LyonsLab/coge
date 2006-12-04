# -*- perl -*-

# t/015_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::SequenceType' ); }

my $object = CoGeX::SequenceType->new ();
isa_ok ($object, 'CoGeX::SequenceType');


