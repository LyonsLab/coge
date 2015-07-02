#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Filter out reads that aren't contained within sequence prior 
#           to cufflinks step to prevent cufflinks errors:
#           "Error (GFaSeqGet): subsequence cannot be larger than 616" and "Error getting subseq for CUFF.26406.1 (1..686)!"
# Author:   mdb
# Created:  6/26/15
#------------------------------------------------------------------------------

use strict;
use warnings;

use CoGe::Accessory::Utils qw(execute);

my $infile = shift;  # SAM input filename
my $outfile = shift; # filtered SAM output filename

my $cmd = "awk 'BEGIN {OFS=\"\t\"} {split(\$6,C,/[0-9]*/); split(\$6,L,/[SMDIN]/); if (C[2]==\"S\") {\$10=substr(\$10,L[1]+1); \$11=substr(\$11,L[1]+1)}; if (C[length(C)]==\"S\") {L1=length(\$10)-L[length(L)-1]; \$10=substr(\$10,1,L1); \$11=substr(\$11,1,L1); }; gsub(/[0-9]*S/,\"\",\$6); print}' $infile > $outfile";

qx{$cmd};
my $cmdStatus = $?;

if ($cmdStatus != 0) {
    print STDERR "error: command failed with rc=$cmdStatus: $cmd";
    exit(-1);
}

exit;
