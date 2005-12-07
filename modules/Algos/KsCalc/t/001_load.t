# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;
use Data::Dumper;
use strict;
BEGIN { use_ok( 'CoGe::Algos::KsCalc' ); }

my $object = CoGe::Algos::KsCalc->new ();
isa_ok ($object, 'CoGe::Algos::KsCalc');

$object->version(6);
$object->name1("at1g01380");
$object->name2("at4g01060");
$object->palign();
my $locs = $object->get_coding_seq_by_name($object->name1);
#print Dumper $locs;
print $object->gaplessP1,"\n";
print $object->gaplessP2,"\n";
print $object->gaplessD1,"\n";
print $object->gaplessD2,"\n";
