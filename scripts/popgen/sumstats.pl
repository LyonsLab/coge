#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:	Calculate pi and print to file
# Author:	Matt Bomhoff
# Created:	9/21/15
#-------------------------------------------------------------------------------

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use CoGe::Algos::PopGen::Diversity;
use CoGe::Algos::PopGen::FileFormats;
use CoGe::Algos::PopGen::Stats qw(nChoose2);
use Tabix;
use Vcf;

my ($VCF_FILE, $GFF_FILE, $BED_FILE, $OUTPUT_FILE, $GENE_NAME, $DEBUG);
GetOptions(
    "vcf=s"   => \$VCF_FILE,
    "gff=s"   => \$GFF_FILE,
    "out=s"   => \$OUTPUT_FILE,
    "gene=s"  => \$GENE_NAME, # optional gene name for testing
    "debug=i" => \$DEBUG
);

my $NUM_INDIVIDUALS = 12 * 2;

# Load gene annotations
print STDERR "Loading GFF\n";
my $pAnnot = loadGFF(file => $GFF_FILE);

# Load VCF
print STDERR "Loading VCF\n";
#my ($pVariants, $pSites) = loadVCF($VCF_FILE);
my $tabix = Tabix->new(-data => $VCF_FILE);

print STDERR "Calculating pi, theta, and Tajima's D\n";
foreach my $type (sort keys %$pAnnot) {
    foreach my $id (sort keys %{$pAnnot->{$type}}) {
        my $mismatch = 0;
        my $compared = 0;
        my $segregating = 0;
        my $sites = 0;
        
        foreach my $feat (sort { $a->{start} <=> $b->{start} } @{$pAnnot->{$type}{$id}}) {
            my $chr   = $feat->{chr};
            my $start = $feat->{start};
            my $end   = $feat->{end};
            unless ($chr && $start && $end) {
                print STDERR "Feature $type $id is missing coordinates\n";
                next;
            }
            
            my $iter = $tabix->query($chr, $start, $end);
            #print STDERR 'query: ', join(' ', $chr, $start, $end), "\n";
            next unless (defined $iter && $iter->get);
            while (my $data = $tabix->read($iter)) {
                my $r = parseVCF($data);
                #print STDERR Dumper $r->{alleles}, "\n";
                if (keys %{$r->{alleles}} > 1) { # variant
                    if (values %{$r->{alleles}} > 2) {
                        print '# ', join("\t", $type, $id, $r->{chr}, $r->{pos}, 'variant', 'skipped (not biallelic)'), "\n" if $DEBUG;
                        next;
                    }
                    my ($count1, $count2) = values %{$r->{alleles}};
                    $mismatch += (($count1 && $count2) ? ($count1 * $count2) : 0);
                    $compared += nChoose2($NUM_INDIVIDUALS);
                    $segregating++;
                    $sites++;
                    print '# ', join("\t", $type, $id, $r->{chr}, $r->{pos}, 'variant', $count1, $count2), "\n" if $DEBUG;
                }
                else { # invariant
                    $compared += nChoose2($NUM_INDIVIDUALS);
                    $sites++;
                    print '# ', join("\t", $type, $id, $r->{chr}, $r->{pos}, 'invariant'), "\n" if $DEBUG;
                }
            }
        }
        
        my ($pi, $theta, $tajD);
        if ($sites) {
            $pi    = $mismatch / $compared;
            $theta = Theta($segregating, $sites, $NUM_INDIVIDUALS);
            $tajD  = TajimasD($pi, $segregating, $NUM_INDIVIDUALS);
        }
        else {
            #$pi = $theta = $tajD = 'N/A';
            next;
        }
        
        #print $sites, $segregating, $mismatch, $compared, "\n";
        print join("\t", $type, $id, $pi, $theta, $tajD), "\n";
    }
}

exit;
#-------------------------------------------------------------------------------

