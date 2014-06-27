# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;
use Data::Dumper;
BEGIN { use_ok( 'CoGe::Algos::Codeml' ); }

my $file = "./t/test.phy";
my $object = CoGe::Algos::Codeml->new (-alignment=>$file,
						 -save_tempfiles=>1);
isa_ok ($object, 'CoGe::Algos::Codeml');
$object->run;
print STDERR Dumper ($object->results());
