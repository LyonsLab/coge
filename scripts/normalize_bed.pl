#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Normalize the score in a BED file to log10 percent (0<=x<=1).
# Usage:    samtools depth file.bam | normalize_bed.pl
#           Note: used to accept output from "samtools mpileup"
# Author:   mdb
# Created:  2/24/14 for qTeller (Cufflinks) integration pipeline
#-------------------------------------------------------------------------------

use strict;

my $filename = shift;
open(my $fh, $filename);

# Determine max score
my $maxScore = 0;
while (<$fh>) {
    my @col = split("\t");
    my $score = normalize($col[4]);
    $maxScore = $score if ($score > $maxScore);
}

# Generate normalized output
seek($fh, 0, 0);
while (<$fh>) {
    chomp;
    my @col = split("\t");
    if ($maxScore > 0) {
        push @col, $col[4]; # copy original score value into column 7 - this is not standard BED format
        $col[4] = normalize($col[4]) / $maxScore;
    }
    print join("\t", @col), "\n";
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
