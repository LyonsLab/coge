# -*- perl -*-

#
# t/008_genomic_sequence_from_feature.t
#

use Test::More tests => 4;

#1
BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

#2
isa_ok ($s, 'CoGeX');
my $feature = $s->resultset('Feature')->find(1);
my $genomeseq = substr($feature->genomic_sequence(), 0, 10);

#3
is( $genomeseq, "AGACAGAATC");

#4
$feature = $s->resultset('Feature')->find(2);
my @seqs = $feature->genomic_sequence();

is( scalar @seqs, 9);
