#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Strip "gi|" and "lcl|" off of input filename.  Used by Cufflinks pipeline.
# Author:   mdb
# Created:  2/10/14
#------------------------------------------------------------------------------
my @filter = ('gi\|', 'lcl\|');
#-------------------------------------------------------------------------------

use strict;

die "Usage:  fasta_reheader.pl <input_file> <output_file>\n" if ($#ARGV+1 < 2);

my $INPUT_FILE = shift @ARGV;
my $OUTPUT_FILE = shift @ARGV;
open(my $inf, $INPUT_FILE) or
    die("Error: cannot open file '$INPUT_FILE'\n");
open(my $outf, ">$OUTPUT_FILE") or
    die("Error: cannot open file '$OUTPUT_FILE'\n");

while (<$inf>) {
    if (/^\>\s*(\S+)/) {
        my $name = $1;
        foreach (@filter) {
            $name =~ s/^$_//;
        }
        print $outf '>', $name, "\n";
    }
    else {
        print $outf $_;
    }
}

close($inf);
close($outf);
exit;
