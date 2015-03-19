#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Filter SNPs from pileup input into VCF output
# Author:   Matt Bomhoff
# Created:  5/13/14
#------------------------------------------------------------------------------

use warnings;
use strict;

use Getopt::Long qw(GetOptions);

our ($min_read_depth, $min_allele_freq, $min_allele_count, $min_base_quality, $quality_scale);

GetOptions(
    "min_read_depth=s"   => \$min_read_depth,
    "min_allele_freq=s"  => \$min_allele_freq,
    "min_allele_count=s" => \$min_allele_count,
    "min_base_quality=s" => \$min_base_quality,
    "quality_scale=s"    => \$quality_scale,
);

# Set defaults
$min_allele_freq = 0.1 unless (defined $min_allele_freq); # min allele frequency
$min_allele_count = 4 unless (defined $min_allele_count); # min high-quality depth of allele
$min_read_depth = 10 unless (defined $min_read_depth);    # overall min depth, regardless of quality
$quality_scale = 32 unless (defined $quality_scale);      # scale for FASTQ encoding
unless (defined $min_base_quality) {
    $min_base_quality = 20; # definition of HQ (High Quality)
    $min_base_quality += $quality_scale + 1; 
}

# Print VCF header
#TODO

my $parsex = qr/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/;
my $line;
while ($line = <STDIN>) {
    chomp $line;
    my ($ref, $pos, $refAllele, $depth, $seq, $qual) = $line =~ $parsex;
    next unless ($depth and $depth >= $min_read_depth);

    my ($hqDepth, $pAlleles) = getHQAlleles($refAllele, \$seq, \$qual);

    # Count alleles that meet minimum criteria
    my %filtered;
    foreach my $b (keys %$pAlleles) {
        next if ($b eq $refAllele); # skip reference allele
        my $freq = $pAlleles->{$b} / $hqDepth;
        if ($freq >= $min_allele_freq and $pAlleles->{$b} >= $min_allele_count) {
            $filtered{$b} = { freq => $freq, depth => $pAlleles->{$b} };
        }
    }

    if (keys %filtered > 1) {
        foreach my $b (keys %filtered) {
            my $info = 'DP='.$filtered{$b}{depth}.';AF='.$filtered{$b}{freq};
            print join("\t", $ref, $pos, '.', $refAllele, $b, '.', '.', $info), "\n";
        }
    }
}

exit;

#-------------------------------------------------------------------------------

# Parse out the alleles from the seq/qual fields in a pileup line.
sub getHQAlleles {
    my $refbase = shift;
    my $pseq = shift;
    my $pqual = shift;

    my %bases;
    my $count = 0;

    my @as = split(//, $$pseq);
    my @aq = unpack("C*", $$pqual);

    my ($i, $j);
    for (($i, $j) = (0, 0);  $i < length($$pseq);  $i++) {
        my $c = $as[$i];
        die "error 1: $i $j $$pseq $$pqual\n" if (not defined $c);
        if ($c eq '>' or $c eq '<') { # reference skip
            $j++;
            next;
        }
        elsif ($c eq '$') { # end of read
            next;
        }
        elsif ($c eq '^') { # start of read followed by encoded quality
            $i++;
            next;
        }
        elsif ($c eq '+' or $c eq '-') { # indel
            $c = $as[$i+1];
            if ($c =~ /[0-9]/) {
                $i++;
                my $c2 = $as[$i+1];
                if ($c2 =~ /[0-9]/) {
                    my $n = int("$c$c2");
                    $i += $n + 1;
                }
                else {
                    $i += $c;
                }
            }
            next;
        }

        my $q = $aq[$j++];
        die "error 2: $i $j $$pseq $$pqual\n" if (not defined $q);
        if ($q >= $min_base_quality and $c ne 'N') {
            $c = $refbase if ($c eq '.' or $c eq ',');
            $c = uc($c);
            next if ($c =~ /[^ACGT]/);
            $bases{$c}++;
            $count++;
        }
    }
    die "error 3: $i $j $$pseq $$pqual" if ($i != @as or $j != @aq);

    return ($count, \%bases);
}
