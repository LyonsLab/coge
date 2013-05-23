#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use Bio::DB::Sam;

my $bamfile = shift or die "Usage: $0 <bamfile>\nOutput a chrom_sizes file.\n";

my $bam = Bio::DB::Sam->new(-bam=>$bamfile) or die "Couldn't open $bamfile";
my @seqids = $bam->seq_ids;
for my $s (sort @seqids) {
    my $len = $bam->length($s);
    print join("\t",$s,$len),"\n";
}
