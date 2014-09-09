#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Convert GFF file to BED format and log normalize the scores to (0, 1).
#           Skip entries of type "exon".
#           Extract "name" and "fpkm" values from "attributes" column.
# Usage:    perl gff_to_bed.pl <input_file> > <output_file>
# Author:   Matt Bomhoff
# Created:  9/9/14 for the McCarthy Lab
#-------------------------------------------------------------------------------
# Parameters:
my $LOG_TRANSFORM = 0; # Set to 1 to log normalize, or 0 to normalize to percent value.
#-------------------------------------------------------------------------------

use strict;

my $filename = shift;
die "Usage:  perl gff_to_bed.pl <input_file> > <output_file>\n" unless (defined $filename);
open(my $fh, $filename) or 
    die("Error: cannot open file '$filename'\n");

# Determine max normalized score
my $lineNum = 0;
my $maxScore = 0;
while (<$fh>) {
    $lineNum++;
    my @col = split("\t");
    my $score = $col[5];
    die "Error: missing score value at line $lineNum\n" if ($score eq '.');
    $score = normalize($score, $LOG_TRANSFORM);
    $maxScore = $score if ($score > $maxScore);
}

# Convert to BED format
seek($fh, 0, 0);
$lineNum = 0;
while (<$fh>) {
    $lineNum++;
    chomp;
    my ($seqid, undef, $type, $start, $end, $score, $strand, undef, $attr) = split("\t");
    next if $type eq 'exon';
    
    my ($name) = $attr =~ /ID=(\w+)/i;
    my ($fpkm) = $attr =~ /fpkm=([\.\d]+)/i;
    
    if ($maxScore > 0) {
        $score = normalize($score, $LOG_TRANSFORM) / $maxScore;
    }
    
    die "Error: invalid parameters at line $lineNum\n" if (!defined($name) || !defined($fpkm));
    print join("\t", $seqid, $start, $end, $name, $score, $strand, $fpkm), "\n";
}

close($fh);
exit;

#-------------------------------------------------------------------------------
sub normalize {
    my ($n, $logFlag) = @_;
    return log10(abs($n)+1) if ($logFlag); # log normalize
    return $n; # just return input value
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
