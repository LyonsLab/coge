#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Convert Homer peaks output to CoGe CSV format
# Usage:    
# Author:   mdb
# Created:  1/28/16
#-------------------------------------------------------------------------------

use strict;

my $filename = shift;
open(my $fh, $filename);

# Generate normalized output
while (<$fh>) {
    next if /^#/;
    chomp;
    my @col = split("\t");
    my $chr = $col[1];
    my $start = $col[2];
    my $end = $col[3];
    my $strand = ($col[4] eq '+' ? 1 : -1);
    my $count = $col[5];
    my $score = $col[7];
    
    print join(",", $chr, $start, $end, $strand, $count, $score), "\n";
}

close($fh);
exit;

#-------------------------------------------------------------------------------
