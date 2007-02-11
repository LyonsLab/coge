# -*- perl -*-

#
# t/007_protein_sequence_from_feature.t
#

use Test::More tests => 3;

#1
BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

#2
isa_ok ($s, 'CoGeX');
my $feature = $s->resultset('Feature')->find(18);
my $proteinseq = substr($feature->protein_sequence(), 0, 10);

#3
is( $proteinseq, "MFNFFFTMKD");
