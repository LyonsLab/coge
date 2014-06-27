#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Convert a SAMTools pileup depth to BED format.
# Usage:    samtools -D file.bam | pileup_to_bed.pl
# Author:   mdb
# Created:  2/24/14
#-------------------------------------------------------------------------------

use strict;

# Simple per position conversion
#while (<>) {
#    my ($chr, $pos, undef, $depth) = split("\t", $_);
#    if ($depth >= $MIN_DEPTH) {
#        if ($FIX_CHR_NAME) {
#            $chr =~ s/^lcl\|//;
#            $chr =~ s/^gi\|//;
#        }
#        print join("\t", $chr, $pos, $pos, '.', $depth, '+'), "\n";
#    }
#}

# Create intervals
my ($start, $stop);
my ($prevChr, $prevDepth);
while (<>) {
    chomp;
    my ($chr, $pos, undef, $depth) = split("\t");

    if (!defined $start || ($pos-$stop > 1) || $depth != $prevDepth || $chr ne $prevChr) {
        # Print interval
        if (defined $start) {
            print_line();
        }

        # Reset interval
        $start = $pos;
        $prevDepth = $depth;
        $prevChr = $chr;
    }

    $stop = $pos;
}

# Print last interval
if (defined $start and $start != $stop) {
    print_line();
}

exit;

#-------------------------------------------------------------------------------
sub print_line {
    $prevChr =~ s/^lcl\|//;
    $prevChr =~ s/^gi\|//;
    print join("\t", $prevChr, $start, $stop+1, '.', $prevDepth, '+'), "\n";
}
