#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Normalize the score in a BED file to log10 percent (0<=x<=1).
# Usage:    samtools -D file.bam | pileup_to_bed.pl
# Author:   mdb
# Created:  2/24/14 for qTeller integration pipeline
#-------------------------------------------------------------------------------

use strict;

my $filename = shift;
open(my $fh, $filename);

# Determine max score
my $maxScore = 0;
while (<$fh>) {
    my @col = split("\t", $_);
    my $score = normalize($col[4]);
    $maxScore = $score if ($score > $maxScore);
}

# Generate normalized output
seek($fh, 0, 0);
while (<$fh>) {
    my @col = split("\t", $_);
    if ($maxScore > 0) {
        $col[4] = normalize($col[4]) / $maxScore;
    }
    print join("\t", @col);
}

close($fh);
exit;

#-------------------------------------------------------------------------------
sub normalize {
    my $n = shift;
    return log10(abs($n)+1);
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
