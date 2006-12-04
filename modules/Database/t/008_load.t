# -*- perl -*-

# t/008_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::DataInformation' ); }

my $object = CoGeX::DataInformation->new ();
isa_ok ($object, 'CoGeX::DataInformation');


