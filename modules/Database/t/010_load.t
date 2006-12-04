# -*- perl -*-

# t/010_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX::DataSource' ); }

my $object = CoGeX::DataSource->new ();
isa_ok ($object, 'CoGeX::DataSource');


