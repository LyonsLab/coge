# -*- perl -*-

# t/026_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::FeatureName' ); }

my $object = CoGeX::FeatureName->new ();
isa_ok ($object, 'CoGeX::FeatureName');


