#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Sanitize the headers in the given FASTA file.  Used by methylation pipeline.
# Author:   mdb
# Created:  12/3/15
#------------------------------------------------------------------------------

use strict;
use CoGe::Core::Genome qw(fix_chromosome_id);

die "Usage:  fix_chromosome_id.pl <input_file> <output_file>\n" if ($#ARGV+1 < 2);

my $INPUT_FILE = shift @ARGV;
my $OUTPUT_FILE = shift @ARGV;
open(my $inf, $INPUT_FILE) or
    die("Error: cannot open file '$INPUT_FILE'\n");
open(my $outf, ">$OUTPUT_FILE") or
    die("Error: cannot open file '$OUTPUT_FILE'\n");

while (<$inf>) {
    if (/^\>\s*(\S+)/) {
        my $name = fix_chromosome_id($1);
        print $outf '>', $name, "\n";
    }
    else {
        print $outf $_;
    }
}

close($inf);
close($outf);
exit;
