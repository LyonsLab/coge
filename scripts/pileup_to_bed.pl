#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Convert a SAMTools pileup depth to BED format.
# Usage:    samtools -D file.bam | pileup_to_bed.pl
# Author:   mdb
# Created:  2/24/14
#-------------------------------------------------------------------------------
my $MIN_DEPTH = 5;
my $FIX_CHR_NAME = 1;
#-------------------------------------------------------------------------------

use strict;

while (<>) {
    my ($chr, $pos, undef, $depth) = split("\t", $_);
    if ($depth >= $MIN_DEPTH) {
        if ($FIX_CHR_NAME) {
            $chr =~ s/^lcl\|//;
            $chr =~ s/^gi\|//;
        }
        print join("\t", $chr, $pos, $pos, '.', $depth, '+'), "\n";
    }
}

exit;

#-------------------------------------------------------------------------------
