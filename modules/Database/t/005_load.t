# -*- perl -*-

# t/005_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::GenomicSequence' ); }

my $object = CoGeX::GenomicSequence->new ();
isa_ok ($object, 'CoGeX::GenomicSequence');


