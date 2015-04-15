#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::Web;
use CoGeX::Result::Genome qw( get_genome_seq );
use Data::Dumper;
use CGI;

my $q = CGI->new;
my $gid = $q->param('gid');
my $chr = $q->param('chr');

# Connect to the database
my ( $db, $user, $conf ) = CoGe::Accessory::Web->init();

my $genome = $db->resultset('Genome')->find($gid);
my $chromosome = $db->resultset('GenomicSequence')->find({genome_id=>$gid,chromosome=>$chr});

print "Content-Type: application/force-download	\n";
print "Content-disposition: attachement; filename=chromosome_" . $gid . "_" . $chr . ".faa\n\n";

print ">chromosome " . $chr . "\n";

# Get sequence from file
my $seq = $genome->get_genomic_sequence(
	chr   => $chr,
	start => 1,
	stop  => $chromosome->sequence_length
);

for (my $i=0; $i<length($seq); $i+=70) {
	print substr($seq, $i, 70);
	print "\n";
}