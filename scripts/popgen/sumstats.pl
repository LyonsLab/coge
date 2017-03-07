#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:	Calculate pi, theta, and Tajima's D and print to STDOUT
# Author:	Matt Bomhoff
# Created:	9/21/15
#-------------------------------------------------------------------------------
# Installing Tabix:
#    git clone git@github.com:samtools/tabix.git
#    cd tabix
#    make
#    sudo cp bgzip tabix /usr/local/bin
#    cd perl
#    perl Makefile.PL lib=/usr/local/lib/perl/5.18.2/
#    sudo make install
#-------------------------------------------------------------------------------

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);
use File::Touch;
use CoGe::Algos::PopGen::Diversity;
use CoGe::Algos::PopGen::FileFormats;
use lib '/opt/apache2/coge/bin/Tabix'; # mdb added 5/4/16 to get Tabix.pm working on PROD
use Tabix;

my ($VCF_FILE, $GFF_FILE, $FASTA_FILE, $CHR, $GENE_NAME, $FEAT_TYPE, 
    $OUTPUT_PATH, $DEBUG, $DEBUGFH);
    
GetOptions(
    "vcf=s"    => \$VCF_FILE,     # required VCF input file
    "gff=s"    => \$GFF_FILE,     # required GFF input file
    "fasta=s"  => \$FASTA_FILE,   # required reference FASTA input file
    "output=s" => \$OUTPUT_PATH,  # optional output path for results files
    "gene=s"   => \$GENE_NAME,    # optional gene name for testing
    "type=s"   => \$FEAT_TYPE,    # optional GFF type for testing
    "chr=s"    => \$CHR,          # optional chromosome name for testing
    "debug"    => \$DEBUG         # optional debug output
);

my $NUM_INDIVIDUALS = 12 * 2;
my $IGNORE_CHROMOSOMES = 1;

# Create output path if specified
if (defined $OUTPUT_PATH) {
    make_path($OUTPUT_PATH) unless (-r $OUTPUT_PATH);
    unless (-r $OUTPUT_PATH) {
        print STDERR "Error: couldn't create output path '$OUTPUT_PATH'\n";
        exit(-1);
    }
}
else {
    $OUTPUT_PATH = './';
}

# Create debug file
if ($DEBUG) {
    my $logfile = catfile($OUTPUT_PATH, 'sumstats.log');
    open($DEBUGFH, ">$logfile") or die "Cannot open log output file '$logfile'\n";
}

# Load GFF
print STDERR "Loading GFF\n";
my $pAnnot = loadGFF(file => $GFF_FILE);

# Load FASTA
print STDERR "Loading FASTA\n";
my $pSeq = loadFASTA($FASTA_FILE);

# Open Tabix indexed VCF
print STDERR "Loading VCF\n";
my $tabix = Tabix->new(-data => $VCF_FILE);

# Create results file
my $resultsfile = catfile($OUTPUT_PATH, 'sumstats.tsv');
open(my $fh, ">$resultsfile") or die "Cannot open results output file '$resultsfile'\n";

# Iterate through GFF entities and calculate summary stats
print STDERR "Calculating pi, theta, and Tajima's D\n";
my %allele_cache;
foreach my $type (sort keys %$pAnnot) {
    next if (defined $FEAT_TYPE and $type ne $FEAT_TYPE);
    next if ($IGNORE_CHROMOSOMES and $type eq 'chromosome');
    
    # Print type header line
    print $fh join("\t", "#$type", 'CHROMOSOME', 'GENE NAME', 'START', 'END', 
        'TOTAL SITES', 'TOTAL SEG. SITES', 'TOTAL PI', 'TOTAL THETA', "TOTAL TAJIMA'S D");
    print $fh "\t", join("\t", 
        '0-FOLD SITES', '0-FOLD SEG. SITES', '0-FOLD PI', '0-FOLD THETA', "0-FOLD TAJIMA'S D", 
        '4-FOLD SITES', '4-FOLD SEG. SITES', '4-FOLD PI', '4-FOLD THETA', "4-FOLD TAJIMA'S D")
        if ($type eq 'cds');
    print $fh "\n";    
    
    foreach my $chr (sort keys %{$pAnnot->{$type}}) {
        next if (defined $CHR and $chr ne $CHR);
    
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
                if ($type eq 'cds') {
                    my $featSeq = substr($pSeq->{$chr}, $start-1, $featLen);
                    debug("$id $chr:$start-$end $featSeq");
                    $codingSeq .= $featSeq;
                }
                
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
                    
                    my $degeneracy;
                    if ($type eq 'cds') {
                        my $codingPos = $pos - $start + $codingOffset;
                        my $codonOffset = $codingPos % 3;
                        my $codonStart = $codingPos - $codonOffset;
                        my $codonSeq = substr($codingSeq, $codonStart, 3);
                        $degeneracy = getDegeneracy($codonSeq, $codonOffset);
                        debug("degeneracy=", defined $degeneracy ? $degeneracy : '?', " codingPos=$codingPos codonSeq=$codonSeq codonOffset=$codonOffset");
                    }
                    
                    my $numAlleles = scalar keys %$alleles;
                    if ($numAlleles < 2) { # invariant
                        debug(join("\t", "$chr:$pos", 'invariant', $type, $id));
                    }
                    else { # variant
                        # Ignore sites that are not biallelic
                        if ($numAlleles > 2) {
                            debug(join("\t", "$chr:$pos", 'variant', $type, $id, 'skipped (not biallelic)'));
                            next;
                        }
                        
                        my ($count1, $count2) = values %$alleles;
                        foreach my $d ('total', $degeneracy) {
                            next unless defined $d;
                            $stats{$d}{mismatch} += (($count1 && $count2) ? ($count1 * $count2) : 0);
                            $stats{$d}{segregating}++;
                            debug(join("\t", "$chr:$pos", 'variant', $type, $id, $count1, $count2));
                        }
                    }
                    
                    foreach my $d ('total', $degeneracy) {
                        next unless defined $d;
                        $stats{$d}{compared} += nChoose2($NUM_INDIVIDUALS);
                        $stats{$d}{sites}++;
                    }
                }
                
                $codingOffset += $featLen;
            }
            
            # Skip if no coverage
            next unless $stats{total}{sites};
            
            # Calculate summary stats
            foreach my $d ('total', 0, 4) {
                my $sites = $stats{$d}{sites};
                next unless $sites;
                my $mismatch = $stats{$d}{mismatch} // 0;
                my $compared = $stats{$d}{compared}; 
                my $segregating = $stats{$d}{segregating} // 0;
                
                $stats{$d}{pi}    = $mismatch / $compared;
                $stats{$d}{theta} = Theta($segregating, $sites, $NUM_INDIVIDUALS);
                $stats{$d}{tajD}  = TajimasD($stats{$d}{pi}, $segregating, $NUM_INDIVIDUALS);
            }
            
            # Print chromosome header line
            print $fh "#$chr\t", length($pSeq->{$chr}), "\n";
            
            # Print result row
            print $fh join("\t", $chr, $id, $featStart, $featEnd);
            my @output = ( $type eq 'cds' ? ('total', 0, 4) : ('total') );
            foreach my $d (@output) {
                print $fh "\t", join("\t", ( $stats{$d}{sites}   // 0,
                                   $stats{$d}{segregating} // 0,
                                   $stats{$d}{pi}    ? sprintf("%.4f", $stats{$d}{pi})    : 0, 
                                   $stats{$d}{theta} ? sprintf("%.4f", $stats{$d}{theta}) : 0, 
                                   $stats{$d}{tajD}  ? sprintf("%.4f", $stats{$d}{tajD})  : 'NaN' ) );
            }
            print $fh "\n";
        }
    }
}

close($DEBUGFH) if $DEBUGFH;

# Create "log.done" file to indicate completion to JEX
my $logdonefile = catfile($OUTPUT_PATH, 'sumstats.done');
touch($logdonefile);

exit;

#-------------------------------------------------------------------------------
sub debug {
    print $DEBUGFH @_, "\n" if $DEBUGFH;
}
