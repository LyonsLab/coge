# -*- perl -*-

# t/015_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'CoGe::Genome::DB::User_group' ); }

my $object = CoGe::Genome::DB::User_group->new ();
isa_ok ($object, 'CoGe::Genome::DB::User_group');


