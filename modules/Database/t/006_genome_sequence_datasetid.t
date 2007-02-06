# -*- perl -*-

# t/006_genome_sequence_datasetid.t
#

use Test::More tests => 6;

BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

#$s->storage->debug(1);

isa_ok ($s, 'CoGeX');
my $dataset = $s->resultset('Dataset')->find(193); # 193 is the dataset_id for
                                                   # chr1 from TIGR's rice,
                                                   # version 4

is( length($dataset->get_genome_sequence( 1, 40, 2239 )), 2200);
is( length($dataset->get_genome_sequence( 1, 40, 12039 )), 12000);
is( length($dataset->get_genome_sequence( 1, 40, 27039 )), 27000);
is( length($dataset->get_genome_sequence( 1, 10, 19 )), 10);
