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

# Determine max normalized score & save mRNA FPKMs
my $lineNum = 0;
my $maxScore = 0;
my %mrnaScores;
while (<$fh>) {
    $lineNum++;
    my @col = split("\t");
    my $type  = $col[2];
    next unless ($type =~ /mrna/i);
    my $score = $col[5];
    my $attr  = $col[8];
    die "Error: missing score value at line $lineNum\n" if ($score eq '.');
    
    # Save mRNA score & fpkm
    my ($id) = $attr =~ /ID=(\w+)/i;
    my ($fpkm) = $attr =~ /fpkm=([\.\d]+)/i;
    die "Error: mRNA missing ID and FPKM attributes at line $lineNum\n" unless (defined($id) && defined($fpkm));
    $mrnaScores{$id} = $fpkm;

    # Max score
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
    next unless ($type eq 'exon');
    die "Error: missing score value at line $lineNum\n" unless (defined $score);
    
    my ($parent) = $attr =~ /Parent=(\w+)/i;
    die "Error: exon is missing PARENT attribute at line $lineNum\n" unless (defined $parent);

    if ($maxScore > 0) {
        $score = normalize($score, $LOG_TRANSFORM) / $maxScore;
    }

    my $fpkm = $mrnaScores{$parent};
    
    print join("\t", $seqid, $start, $end, $parent, $score, $strand, $fpkm), "\n";
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
