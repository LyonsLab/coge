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
use Tabix;
use Vcf;

my ($VCF_FILE, $GFF_FILE, $FASTA_FILE, $CHR, $GENE_NAME, $FEAT_TYPE, $DEBUG_FILE, $DEBUGFH);
GetOptions(
    "vcf=s"   => \$VCF_FILE,    # required VCF input file
    "gff=s"   => \$GFF_FILE,    # required GFF input file
    "fasta=s" => \$FASTA_FILE,  # required reference FASTA input file
    "gene=s"  => \$GENE_NAME,   # optional gene name for testing
    "type=s"  => \$FEAT_TYPE,   # optional GFF type for testing
    "chr=s"   => \$CHR,         # optional chromosome name for testing
    "debug=s" => \$DEBUG_FILE   # optional debug output filename
);

my $NUM_INDIVIDUALS = 12 * 2;

open($DEBUGFH, ">$DEBUG_FILE") or die "Cannot open debug output file\n" if $DEBUG_FILE;

# Load GFF
print STDERR "Loading GFF\n";
my $pAnnot = loadGFF(file => $GFF_FILE);

# Load FASTA
print STDERR "Loading FASTA\n";
my $pSeq = loadFASTA($FASTA_FILE);

# Load VCF
print STDERR "Loading VCF\n";
#my ($pVariants, $pSites) = loadVCF($VCF_FILE);
my $tabix = Tabix->new(-data => $VCF_FILE);

# Iterate through GFF entities and calculate summary stats
print STDERR "Calculating pi, theta, and Tajima's D\n";
my %allele_cache;
foreach my $type (sort keys %$pAnnot) {
    next if (defined $FEAT_TYPE and $type ne $FEAT_TYPE);
    
    foreach my $chr (sort keys %{$pAnnot->{$type}}) {
        next if (defined $CHR and $chr ne $CHR);
    
        print join("\t", "#$type $chr", 'GENE NAME', 'START', 'END', 'TOTAL SITES', 'SEG. SITES', 'COMPARISONS', 'MISMATCHES', 'PI', 'THETA', "TAJIMA'S D"), "\n";
        foreach my $id (sort keys %{$pAnnot->{$type}{$chr}}) {
            next if (defined $GENE_NAME and $id ne $GENE_NAME);
            
            my %stats;
            #my ($sites, $segregating, $compared, $mismatch) = (0, 0, 0, 0);
            my ($codingOffset, $codingSeq, $featStart, $featEnd, $featChr);
            $codingOffset = 0;
            
            foreach my $feat (sort { $a->{start} <=> $b->{start} } @{$pAnnot->{$type}{$chr}{$id}}) {
                my $start = $feat->{start};
                my $end   = $feat->{end};
                unless ($chr && $start && $end) {
                    print STDERR "Feature $type $id is missing coordinates\n";
                    next;
                }
                
                $featStart = $start unless defined $featStart;
                $featEnd = $end if ((not defined $featEnd) || $end > $featEnd);
                my $featLen = $end - $start + 1;
                my $featSeq = substr($pSeq->{$chr}, $start-1, $featLen);
                debug("$id $chr:$start-$end $featSeq");
                $codingSeq .= $featSeq;
                
                # Retrieve range of entries from indexed VCF file
                my $iter = $tabix->query($chr, $start, $end);
                #print STDERR 'query: ', join(' ', $chr, $start, $end), "\n";
                next unless (defined $iter && $iter->get);
                
                # Count and classify sites
                while (my $data = $tabix->read($iter)) {
                    my $r = parseVCF($data);
                    my $pos = $r->{pos};
                    
                    unless ($allele_cache{$chr} && $allele_cache{$chr}{$pos}) {
                        $allele_cache{$chr}{$pos} = parseGenotypes($r);
                    }
                    my $alleles = $allele_cache{$chr}{$pos};
                    #print STDERR Dumper $alleles, "\n";
                    
                    my $codingPos = $pos - $start + $codingOffset;
                    my $codonOffset = $codingPos % 3;
                    my $codonStart = $codingPos - $codonOffset;
                    my $codonSeq = substr($codingSeq, $codonStart, 3);
                    my $degeneracy = getDegeneracy($codonSeq, $codonOffset);
                    debug("degeneracy=$degeneracy codingPos=$codingPos codonSeq=$codonSeq codonOffset=$codonOffset");
                    
                    my $numAlleles = scalar keys %$alleles;
                    if ($numAlleles < 2) { # invariant
                        debug('# ', join("\t", $type, $id, $chr, $pos, 'invariant'));
                    }
                    else { # variant
                        if ($numAlleles > 2) {
                            debug('# ', join("\t", $type, $id, $chr, $pos, 'variant', 'skipped (not biallelic)'));
                            next;
                        }
                        
                        my ($count1, $count2) = values %$alleles;
                        foreach my $d ('all', $degeneracy) {
                            $stats{$d}{mismatch} += (($count1 && $count2) ? ($count1 * $count2) : 0);
                            $stats{$d}{segregating}++;
                            debug('# ', join("\t", $type, $id, $chr, $pos, 'variant', $count1, $count2));
                        }
                    }
                    
                    foreach my $d ('all', $degeneracy) {
                        $stats{$d}{compared} += nChoose2($NUM_INDIVIDUALS);
                        $stats{$d}{sites}++;
                    }
                }
                
                $codingOffset += $featLen;
            }
            
            # Skip if no coverage
            next unless $stats{all}{sites};
            
            # Calculate summary stats
            foreach my $degeneracy ('all', 0, 4) {
                next unless $stats{$degeneracy}{sites};
                $stats{$degeneracy}{pi}    = $stats{$degeneracy}{mismatch} / $stats{$degeneracy}{compared};
                $stats{$degeneracy}{theta} = Theta($stats{$degeneracy}{segregating}, $stats{$degeneracy}{sites}, $NUM_INDIVIDUALS);
                $stats{$degeneracy}{tajD}  = TajimasD($stats{$degeneracy}{pi}, $stats{$degeneracy}{segregating}, $NUM_INDIVIDUALS);
            }
            
            # Print to output
            #print join("\t", $id, $featStart, $featEnd, $sites, $segregating, $compared, $mismatch, $pi, $theta, $tajD), "\n";
            print join("\t", $id, $featStart, $featEnd);
            foreach my $degeneracy ('all', 0, 4) {
                print join("\t", ($stats{$degeneracy}{sites} // 0, 
                                 $stats{$degeneracy}{segregating} // 0, 
                                 $stats{$degeneracy}{compared} // 0,  
                                 $stats{$degeneracy}{mismatch} // 0, 
                                 $stats{$degeneracy}{pi} // '', 
                                 $stats{$degeneracy}{theta} // '', 
                                 $stats{$degeneracy}{tajD} // '') );
            }
            print "\n";
        }
    }
}

close($DEBUGFH) if $DEBUGFH;

exit;
#-------------------------------------------------------------------------------
sub debug {
    print $DEBUGFH @_, "\n" if $DEBUGFH;
}
