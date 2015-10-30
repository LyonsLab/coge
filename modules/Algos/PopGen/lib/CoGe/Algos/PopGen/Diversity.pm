#-------------------------------------------------------------------------------
# Purpose:	Functions for detecting/quantifying diversity and divergence.
# Author:	Matt Bomhoff
# Created:	7/7/11, imported into CoGe 9/21/15
#-------------------------------------------------------------------------------
package CoGe::Algos::PopGen::Diversity;

use warnings;
use strict;
use base 'Exporter';
use List::Util qw(first);

our $DEBUG = 0;

our @EXPORT = qw($DEBUG isStopCodon getAA getDegeneracy getSynonymity 
				 getBasesAtPos getSites subSampleSites maskSequence 
				 getMutations getMutationsForSubsampledSites getPolyWithSubsampling  
				 JukesCantor Pi PiForSites PiForPolySites PiForSitesWithSubsampling
				 PiForBiallelic Theta PairwiseMismatches
				 MeanPairwiseDiff DxyForSites  DxyForPolySites Fst FstForSites KaKsForSites 
				 TajimasD FuLiD harmonic nChoose2
				 %SPECIES @SUB_LINES @ALL_LINES);

our %SPECIES = (
	'CAS' => [ 'CIM', 'CKN', 'CKS', 'CTP', 'DKN', 'MDG', 'MPR', 'sanger_CAST' ],
	'CAS6' => [ 'CIM', 'CKN', 'CKS', 'DKN', 'MDG', 'sanger_CAST' ],
	'DOM' => [ 'BIK', 'BZO', 'DCP', 'DJO', 'DMZ', 'LEWES', 'WLA', 'sanger_WSB' ], 
	'MUS' => [ 'BID', 'CZCH2', 'MBK', 'MBT', 'MCZ', 'MDH', 'MPB', 'merged_PWK' ],
	'MUS7' => [ 'CZCH2', 'MBK', 'MBT', 'MCZ', 'MDH', 'MPB', 'merged_PWK' ],
	'CAROLI' => [ 'CAROLI' ],
	'SPRET' => [ 'merged_SPRET' ],
	
	'CAS6_DOM_MUS7' => [ 'CIM', 'CKN', 'CKS', 'DKN', 'MDG', 'sanger_CAST', 
						 'BIK', 'BZO', 'DCP', 'DJO', 'DMZ', 'LEWES', 'WLA', 'sanger_WSB', 
						 'CZCH2', 'MBK', 'MBT', 'MCZ', 'MDH', 'MPB', 'merged_PWK' ] 
);

our @SUB_LINES = (@{$SPECIES{'CAS6'}}, @{$SPECIES{'DOM'}}, @{$SPECIES{'MUS7'}});
our @ALL_LINES = (@SUB_LINES, 'merged_SPRET', 'CAROLI');

my %codonTable = (
	'AAT' => 'N', 'AAC' => 'N', # Asparagine
	'AAA' => 'K', 'AAG' => 'K', # Lysine
	'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T', 'ACN' => 'T', 	# Threonine
	'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', # Isoleucine
	'ATG' => 'M', 				# Methionine
	'CAT' => 'H', 'CAC' => 'H', # Histidine
	'CAA' => 'Q', 'CAG' => 'Q', # Glutamine
	'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'CGN' => 'R', 'AGA' => 'R', 'AGG' => 'R', # Arginine
	'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P', 'CCN' => 'P', # Proline
	'GAT' => 'D', 'GAC' => 'D', # Aspartic Acid
	'GAA' => 'E', 'GAG' => 'E', # Glutamic Acid
	'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'GCN' => 'A', # Alanine
	'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'GGN' => 'G', # Glycine
	'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'GTN' => 'V', # Valine
	'AGT' => 'S', 'AGC' => 'S', 'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'TCN' => 'S', # Serine
	'TAT' => 'Y', 'TAC' => 'Y', # Tyrosine
	'TAA' => '_', 'TAG' => '_', 'TGA' => '_', # Stop
	'TTT' => 'F', 'TTC' => 'F', # Phenylalanine
	'TTA' => 'L', 'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'CTN' => 'L', # Leucine
	'TGT' => 'C', 'TGC' => 'C', # Cysteine
	'TGG' => 'W' 				# Tryptophan
);

my %foldTable = ( # mdb updated 8/9/11
	0 => {
		'NCN' => 0, 'NNC' => 0, 'NNT' => 0,
		'NAC' => 0, 'NAT' => 0,
		'NCA' => 0, 'NCC' => 0, 'NCG' => 0, 'NCT' => 0,
		'AAA' => 0, 'AAC' => 0, 'AAG' => 0, 'AAT' => 0, 'AAN' => 0, # Lysine/Asparagine
	    'ACA' => 0, 'ACC' => 0, 'ACG' => 0, 'ACT' => 0, 'ACN' => 0, # Threonine
	    'AGA' => 2, 'AGC' => 0, 'AGG' => 2, 'AGT' => 0, 			# Arginine/Serine
	    'ATA' => 2, 'ATC' => 0, 'ATG' => 2, 'ATT' => 0, 			# Isoleucine/Methionine
	    'CAA' => 0, 'CAC' => 0, 'CAG' => 0, 'CAT' => 0, 'CAN' => 0, # Glutamine/Histidine
	    'CCA' => 0, 'CCC' => 0, 'CCG' => 0, 'CCT' => 0, 'CCN' => 0, # Proline
	    'CGA' => 2, 'CGC' => 0, 'CGG' => 2, 'CGT' => 0, 			# Arginine
	    'CTA' => 2, 'CTC' => 0, 'CTG' => 2, 'CTT' => 0, 			# Leucine
	    'GAA' => 0, 'GAC' => 0, 'GAG' => 0, 'GAT' => 0, 'GAN' => 0, # Glutamic Acid/Aspartic Acid
	    'GCA' => 0, 'GCC' => 0, 'GCG' => 0, 'GCT' => 0, 'GCN' => 0, # Alanine
	    'GTA' => 2, 'GTC' => 0, 'GTG' => 2, 'GTT' => 0, 			# Valine
	    'GGA' => 2, 'GGC' => 0, 'GGG' => 2, 'GGT' => 0, 			# Glycine
	    'TAC' => 0, 'TAT' => 0,  									# Tyrosine
	    'TCA' => 0, 'TCC' => 0, 'TCG' => 0, 'TCT' => 0, 'TCN' => 0, # Serine
	    'TGC' => 0, 'TGG' => 2, 'TGT' => 0, 						# Cysteine/Tryptophan
	    'TTA' => 2, 'TTC' => 0, 'TTG' => 2, 'TTT' => 0, 			# Leucine/Phenylalanine
	    'TAA' => 0, 'TAG' => 0, 'TGA' => 2, 						# Stop
	    'TAN' => 0, 'NAA' => 0, 'NAG' => 0, 'NAN' => 0, 
	},
	# position 1 is always 0-fold
	2 => {
		'NCN' => 4,
		'NAA' => 2, 'NAC' => 2, 'NAG' => 2, 'NAT' => 2, 'NAN' => 2,
		'NCA' => 4, 'NCC' => 4, 'NCG' => 4, 'NCT' => 4,
		'AAA' => 2, 'AAC' => 2, 'AAG' => 2, 'AAT' => 2, 'AAN' => 2, # Lysine/Asparagine
		'ACA' => 4, 'ACC' => 4, 'ACG' => 4, 'ACT' => 4, 'ACN' => 4, # Threonine
	    'AGA' => 2, 'AGC' => 2, 'AGG' => 2, 'AGT' => 2, 'AGN' => 2, # Arginine/Serine
	    'ATA' => 2, 'ATC' => 2, 'ATG' => 2, 'ATT' => 2, 'ATN' => 2,	# Isoleucine/Methionine (actually 3-fold, using standard of 2)
	    'CAA' => 2, 'CAG' => 2, 'CAC' => 2, 'CAT' => 2, 'CAN' => 2, # Histidine/Glutamine
	    'CCA' => 4, 'CCC' => 4, 'CCG' => 4, 'CCT' => 4, 'CCN' => 4, # Proline
	   	'CGA' => 4, 'CGC' => 4, 'CGG' => 4, 'CGT' => 4, 'CGN' => 4, # Arginine
	   	'CTA' => 4, 'CTC' => 4, 'CTG' => 4, 'CTT' => 4, 'CTN' => 4, # Leucine 
		'GAA' => 2, 'GAG' => 2, 'GAC' => 2, 'GAT' => 2, 'GAN' => 2, # Aspartic Acid/Glutamic Acid
	    'GCA' => 4, 'GCC' => 4, 'GCG' => 4, 'GCT' => 4, 'GCN' => 4, # Alanine
	    'GGA' => 4, 'GGC' => 4, 'GGG' => 4, 'GGT' => 4, 'GGN' => 4, # Glycine
	    'GTA' => 4, 'GTC' => 4, 'GTG' => 4, 'GTT' => 4, 'GTN' => 4, # Valine
	    'TCA' => 4, 'TCC' => 4, 'TCG' => 4, 'TCT' => 4, 'TCN' => 4, # Serine
	    'TTA' => 2, 'TTG' => 2, 'TTC' => 2, 'TTT' => 2, 'TTN' => 2, # Leucine/Phenylalanine
	   	'TGC' => 2, 'TGG' => 2, 'TGT' => 2,  						# Cysteine/Tryptophan
	    'TAC' => 2, 'TAT' => 2, 									# Tyrosine
	    'TAA' => 2, 'TAG' => 2, 'TGA' => 2, 'TAN' => 2, 'TGN' => 2,	# Stop
	}
);

my %synTable = ( # Polymorphorama values
    'AGC' => (1/3), 'AGT' => (1/3), 'TCA' => 1, 'TCC' => 1, 'TCG' => 1, 'TCT' => 1, 'TCN' => 1, # Serine
    'TTC' => (1/3), 'TTT' => (1/3), # Phenylalanine
    'TTA' => (2/3), 'TTG' => (2/3), 'CTA' => (4/3), 'CTC' => 1, 'CTG' => (4/3), 'CTT' => 1, # Leucine
    'TAC' => 1, 'TAT' => 1, # Tyrosine
    #'TAA' => '_', 'TAG' => '_', 'TGA' => '_', # Stop
    'TGC' => (1/2), 'TGT' => (1/2), # Cysteine
    'TGG' => 0, # Tryptophan
    'CCA' => 1, 'CCC' => 1, 'CCG' => 1, 'CCT' => 1, 'CCN' => 1,  # Proline
    'CAC' => (1/3), 'CAT' => (1/3), # Histidine
    'CAA' => (1/3), 'CAG' => (1/3), # Glutamine
    'AGA' => (5/6), 'AGG' => (2/3), 'CGA' => (4/3), 'CGC' => 1, 'CGG' => (4/3), 'CGT' => 1, 'CGN' => 1, # Arginine
    'ATA' => (2/3), 'ATC' => (2/3), 'ATT' => (2/3), # Isoleucine
    'ATG' => 0, # Methionine
    'ACA' => 1, 'ACC' => 1, 'ACG' => 1, 'ACT' => 1, 'ACN' => 1, # Threonine
    'AAC' => (1/3), 'AAT' => (1/3), # Asparagine
    'AAA' => (1/3), 'AAG' => (1/3), # Lysine
    'GTA' => 1, 'GTC' => 1, 'GTG' => 1, 'GTT' => 1, 'GTN' => 1, # Valine
    'GCA' => 1, 'GCC' => 1, 'GCG' => 1, 'GCT' => 1, 'GCN' => 1, # Alanine
    'GAC' => (1/3), 'GAT' => (1/3), # Aspartic Acid
    'GAA' => (1/3), 'GAG' => (1/3), # Glutamic Acid
    'GGA' => 1, 'GGC' => 1, 'GGG' => 1, 'GGT' => 1, 'GGN' => 1, # Glycine
);

my %synTableByPos = ( # Reviewed by Matt & Megan 9/29/11
	0 => {
	    'GCA' => 0, 'GCC' => 0, 'GCG' => 0, 'GCT' => 0, 'GCN' => 0, # Alanine
	    'TGC' => 0, 'TGT' => 0, # Cysteine
	    'GAC' => 0, 'GAT' => 0, # Aspartic Acid
	    'GAA' => 0, 'GAG' => 0, # Glutamic Acid
	    'TTC' => 0, 'TTT' => 0, # Phenylalanine
	    'GGA' => 0, 'GGC' => 0, 'GGG' => 0, 'GGT' => 0, 'GGN' => 0, # Glycine
		'CAC' => 0, 'CAT' => 0, # Histidine
	    'ATA' => 0, 'ATC' => 0, 'ATT' => 0, # Isoleucine
	    'AAA' => 0, 'AAG' => 0, # Lysine
	    'TTA' => (1/6), 'TTG' => (1/6), 'CTA' => (1/6), 'CTC' => 0, 'CTG' => (1/6), 'CTT' => 0, # Leucine 
	    'AAC' => 0, 'AAT' => 0, # Asparagine
	    'CCA' => 0, 'CCC' => 0, 'CCG' => 0, 'CCT' => 0, 'CCN' => 0, # Proline
	    'CAA' => 0, 'CAG' => 0, # Glutamine
	    'AGA' => (1/6), 'AGG' => (1/6), 'CGA' => (1/6), 'CGC' => 0, 'CGG' => (1/6), 'CGT' => 0, # Arginine
		'AGC' => 0, 'AGT' => 0, 'TCA' => 0, 'TCC' => 0, 'TCG' => 0, 'TCT' => 0, # Serine
	    'ACA' => 0, 'ACC' => 0, 'ACG' => 0, 'ACT' => 0, 'ACN' => 0, # Threonine
	    'GTA' => 0, 'GTC' => 0, 'GTG' => 0, 'GTT' => 0, 'GTN' => 0, # Valine
	    'TAC' => 0, 'TAT' => 0, # Tyrosine
	    'ATG' => 0, # Methionine
	    'TGG' => 0, # Tryptophan
	    #'TGA' => '_', 'TAA' => '_', 'TAG' => '_', # Stop
	},
	2 => {
	    'GCA' => 1, 'GCC' => 1, 'GCG' => 1, 'GCT' => 1, 'GCN' => 1, # Alanine
	    'TGC' => (1/6), 'TGT' => (1/6), # Cysteine
	    'GAC' => (1/3), 'GAT' => (1/3), # Aspartic Acid
	    'GAA' => (1/3), 'GAG' => (1/3), # Glutamic Acid
	    'TTC' => (1/3), 'TTT' => (1/3), # Phenylalanine
	    'GGA' => 1, 'GGC' => 1, 'GGG' => 1, 'GGT' => 1, 'GGN' => 1, # Glycine
		'CAC' => (1/3), 'CAT' => (1/3), # Histidine
	    'ATA' => (1/2), 'ATC' => (1/2), 'ATT' => (1/2), # Isoleucine
	    'AAA' => (1/3), 'AAG' => (1/3), # Lysine
	    'TTA' => (1/3), 'TTG' => (1/3), 'CTA' => 1, 'CTC' => 1, 'CTG' => 1, 'CTT' => 1, 'CTN' => 1, # Leucine 
	    'AAC' => (1/3), 'AAT' => (1/3), # Asparagine
	    'CCA' => 1, 'CCC' => 1, 'CCG' => 1, 'CCT' => 1, 'CCN' => 1, # Proline
	    'CAA' => (1/3), 'CAG' => (1/3), # Glutamine
	    'AGA' => (1/3), 'AGG' => (1/3), 'CGA' => 1, 'CGC' => 1, 'CGG' => 1, 'CGT' => 1, 'CGN' => 1, # Arginine
		'AGC' => (1/3), 'AGT' => (1/3), 'TCA' => 1, 'TCC' => 1, 'TCG' => 1, 'TCT' => 1, 'TCN' => 1, # Serine
	    'ACA' => 1, 'ACC' => 1, 'ACG' => 1, 'ACT' => 1, 'ACN' => 1, # Threonine
	    'GTA' => 1, 'GTC' => 1, 'GTG' => 1, 'GTT' => 1, 'GTN' => 1, # Valine
	    'TAC' => (1/3), 'TAT' => (1/3), # Tyrosine
	    'ATG' => (1/2), # Methionine
	    'TGG' => (1/6), # Tryptophan
	    #'TGA' => '_', 'TAA' => '_', 'TAG' => '_', # Stop
	}
);

# Significance tables for Fu & Li's D and D*
my %D_significance = (
	5 => [ -1.96, -1.77, -1.57, 1.63, 1.83, 2.02 ],
	6 => [ -2.1,  -1.88, -1.69, 1.53, 1.71, 1.88 ],
	7 => [ -2.15, -1.9,  -1.67, 1.48, 1.66, 1.8  ],
	8 => [ -2.27, -1.97, -1.75, 1.45, 1.61, 1.76 ]
);
my %DStar_significance = (
	5 => [ -1.26, -1.23, -1.2,  1.57, 1.68, 1.77 ],
	6 => [ -1.54, -1.49, -1.43, 1.46, 1.55, 1.62 ],
	7 => [ -1.75, -1.67, -1.57, 1.37, 1.46, 1.56 ],
	8 => [ -1.93, -1.82, -1.67, 1.34, 1.43, 1.51 ]
);
#-------------------------------------------------------------------------------

# Converts input codon sequence to an amino acid character
sub getAA {
	my $codon = shift;
	return $codonTable{$codon};
}

sub isStopCodon {
	my $codon = shift;
	return ($codon eq 'TAA' or $codon eq 'TAG' or $codon eq 'TGA');
}

# Returns the degeneracy of a particular site given the codon sequence
sub getDegeneracy {
	my $codon = shift; # string of 3 base characters
	my $pos = shift; # offset into codon (0, 1, 2)
	die if ($pos < 0 or $pos > 2);
	return 0 if ($pos == 1);
	return $foldTable{$pos}{$codon};
}

# Returns the synonymity of the given the codon sequence
#sub getSynonymity {
#	my $codon = shift;
#	return $synTable{$codon};
#}

# Returns the synonymity of a particular site given the codon sequence
sub getSynonymity {
	my $codon = shift; # string of 3 base characters
	my $pos = shift; # offset into codon (0, 1, 2)
	die if ($pos < 0 or $pos > 2);
	return 0 if ($pos == 1);
	return $synTableByPos{$pos}{$codon};
}

sub getBasesAtPos {
	my $pSeq = shift;
	my $pSpecies = shift;
	my $pos = shift;
	
	my $numBases = 0;
	my %bases;
	foreach my $s (@$pSpecies) {
		my $b = substr($pSeq->{$s}, $pos, 1);
		if ($b ne 'N') {
			$bases{$b}++;
			$numBases++;	
		}
	}

	return ($numBases, \%bases);
}

# Exclude:
#    - stop codons
#    - sites with ambiguous degeneracy
#    - sites with >2 alleles across all lines and outgroups
#    - sites with missing coverage > X within subspecies
#    - sites missing outgroup (optional)
#    - singletons (optional)
sub getSites {
	my $pSeq = shift; 				# ref to sequences hashed by name
	my $pSpecies = shift;			# ref to array of species names
	my $outgroup = shift;			# optional: name of outgroup
	my $pFlags = shift;				# optional: hash of flags

	die if (keys %$pSeq == 0);
	die if (not defined $pSpecies or @$pSpecies == 0);
	
	$pFlags = () if (not defined $pFlags);
	$pFlags->{excludeSingletons} = 0 if (not defined $pFlags->{excludeSingletons});
	$pFlags->{requireDegen} 	 = 1 if (not defined $pFlags->{requireDegen});
	$pFlags->{requireSyn} 		 = 0 if (not defined $pFlags->{requireSyn});
	$pFlags->{minBases} 		 = @$pSpecies if (not defined $pFlags->{minBases});
	
	my $len = length( first {defined($_)} values %$pSeq );
	
	my (%sites, %poly, %fixed);#, %freq, %syn, %nonsyn); # sites by fold
	for (my $pos = 0;  $pos < $len;  $pos += 3) { # for each codon
		# Get codon sequences at this position
		my %codons;
		foreach my $name (@$pSpecies, $outgroup) {
			next if (not defined $name); # if no outgroup specified
			die "Sequence for '$name' not found\n" if (not defined $pSeq->{$name});
			my $codon = substr($pSeq->{$name}, $pos, 3);
			$codon .= 'N' x (3-(length $codon)) if (length $codon < 3); # pad to codon boundary
			if (isStopCodon($codon)) {
				print "stop: pos=$pos codon=$codon\n" if ($DEBUG);
				goto NEXT_CODON;
			}
			$codons{$name} = $codon;
		}
		
		for my $i (0, 1, 2) { # for each site within codon
			my $pos2 = $pos + $i;
			
			# Determine polymorphism within subspecies
			my %bases;
			my $numBases = 0;
			foreach my $name (@$pSpecies) {
				my $b = substr($codons{$name}, $i, 1);
				if ($b ne 'N') {
					$numBases++;
					$bases{$b}++; 
				}
			}
			
			# Exclude if too many missing bases - default is no missing allowed
			if ($numBases < $pFlags->{minBases}) {
				print "missing: pos=$pos2 count=$numBases/$pFlags->{minBases}\n" if ($DEBUG);
				goto NEXT_SITE;
			}
			
			# Exclude if missing outgroup (when specified)
			my $outB;
			if (defined $outgroup) {
				die "Sequence for '$outgroup' not found\n" if (not defined $pSeq->{$outgroup});
				$outB = substr($codons{$outgroup}, $i, 1);
				if ($outB eq 'N') {
					print "missing: pos=$pos2 outgroup\n" if ($DEBUG);
					goto NEXT_SITE;
				}
			}
			
			# Exclude singleton polymorphisms (optional)
			if ($pFlags->{excludeSingletons}) {
				my $minorCount = first {defined($_)} sort {$a<=>$b} values %bases;
				if ($minorCount == 1) {
					print "singleton: pos=$pos2\n" if ($DEBUG);
					goto NEXT_SITE;
				}
			}
			
			# Count alleles among all lines, exclude sites with more than 2
			{	my %allBases;
				foreach my $name (@ALL_LINES) {
					my $b = substr($pSeq->{$name}, $pos2, 1);
					$allBases{$b}++ if ($b ne 'N');
				}
				if (keys %allBases > 2) {
					print "alleles: pos=$pos2 codon=$pos bases:" . join(',', keys %allBases) . "\n" if ($DEBUG);
					goto NEXT_SITE;
				}
			}
			
			# Determine site degeneracy, skip site if ambiguous
			my ($fold, $synonymity);
			foreach my $name (@$pSpecies, $outgroup) {
				next if (not defined $name); # if no outgroup specified
			
				my $f = getDegeneracy($codons{$name}, $i);
				$fold = $f if (not defined $fold);
				if ($pFlags->{requireDegen} and (not defined $f or $fold != $f)) {
					print "ambiguous degeneracy: pos=$pos2 codon=$pos $codons{$name}" . (defined $f and defined $fold ? " $fold != $f" : "") . "\n" if ($DEBUG);
					goto NEXT_SITE;
				}
				
#				my $s = getSynonymity($codons{$name}, $i);
#				$synonymity = $s if (not defined $synonymity);
#				if ($pFlags->{requireSyn} and defined $s and $s != $synonymity) {
#					print "ambiguous synonymity: pos=$pos2 codon=$pos $codons{$name} $s != $synonymity\n" if ($DEBUG);
#					goto NEXT_SITE;
#				}
			}
			die "pos=$pos2 codon=$pos" if (not defined $fold);

			# Count amino acids
#			my %aa;
#			foreach my $name (@$pSpecies, $outgroup) {
#				next if (not defined $name); # if no outgroup specified
#				my $a = getAA($codons{$name});
#				$aa{$a}++ if (defined $a);
#			}

			# Record site
			$sites{'all'}{$pos2}++;
			$sites{$fold}{$pos2}++;
#			$syn{$pos2} = $synonymity;
#			$nonsyn{$pos2}++ if (keys %aa > 1);
			
			# Determine polymorphism / fixed difference
			if (keys %bases > 1) { # Polymorphic within subspecies
				# Record polymorphism
				$poly{'all'}{$pos2}++;
				$poly{$fold}{$pos2}++;
				print "poly: pos=$pos2 fold=$fold bases: " . join(',', keys %bases) . "\n" if ($DEBUG);
				
				# Record site frequency - FIXME: this is broken if missing data allowed
#				my $minorCount = first {defined($_)} sort {$a<=>$b} values %bases;
#				$freq{$fold}{$minorCount}++;
#				print "SFS: pos=$pos2 fold=$fold count=$minorCount\n" if ($DEBUG);
			}
			elsif (keys %bases == 1) { # Not polymorphic within subspecies
				# Record site frequency (zero)
#				$freq{$fold}{0}++;
				
				# Fixed difference
				if (defined $outgroup and not defined $bases{$outB}) {
					$fixed{'all'}{$pos2}++;
					$fixed{$fold}{$pos2}++;
					print "fixed: pos=$pos2 fold=$fold bases: " . join(',', keys %bases) . " - $outB\n" if ($DEBUG);
				}
			}
			NEXT_SITE: # skip site
		}
		NEXT_CODON: # skip codon
	}
	
	return (\%sites, \%poly, \%fixed);#, \%freq, \%syn, \%nonsyn);
}

sub subSampleSites {
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	my $numSamples = shift;
	my $outgroup = shift; # optional
	my $excludeSingletons = shift; # optional boolean flag
	my (%sites, %poly, %fixed);
	
	foreach my $pos (keys %{$pSites->{'all'}}) {
		# Get subspecies bases
		my @bases;
		foreach my $s (@$pSpecies) {
			my $b = substr($pSeq->{$s}, $pos, 1);
			push @bases, $b if ($b ne 'N');
		}
		next if (@bases < $numSamples);
		
		# Get outgroup base
		my $outB;
		if (defined $outgroup) {
			$outB = substr($pSeq->{$outgroup}, $pos, 1);
			next if ($outB eq 'N');
		}
		
		# Randomly sample bases from subspecies
		while (@bases > $numSamples) {
			my $i = int(rand(@bases));
			splice @bases, $i, 1;
		}
#		@bases = shuffle(@bases);
#		splice(@bases, 0, $numSamples);
		
		# Count alleles
		my %counts;
		foreach my $b (@bases) {
			$counts{$b}++;
		}
		
		# Exclude singletons
		if ($excludeSingletons) {
			my $minorCount = first {defined($_)} sort {$a<=>$b} values %counts;
			next if ($minorCount == 1);	
		}
		
		# Record types of sites by degeneracy
		foreach my $fold ('all', 0, 4) {
			if (defined $pSites->{$fold}{$pos}) {
				$sites{$fold}{$pos} = \%counts;
				if (keys %counts > 1) { # is site still polymorphic?
					$poly{$fold}{$pos} = \%counts;
				}
				elsif (defined $outgroup and not defined $counts{$outB}) { # fixed difference
					$fixed{$fold}{$pos} = \%counts;
				}
			}
		}
	}
	
	return (\%sites, \%poly, \%fixed);
}

sub maskSequence {
	my $pSeq = shift;
	my $pSites = shift;
	
	my $len = length( first {defined($_)} values %$pSeq );
	
	foreach my $s (keys %$pSeq) {
		for (my $pos = 0;  $pos < $len;  $pos++) {
			if (not defined $pSites->{'all'}{$pos}) {
				substr($pSeq->{$s}, $pos, 1) = 'N';
			}
		}
	}
}

sub getMutations {
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	my $outgroup = shift;
	
	my $N = 0;
	my $Ns = 0;
	my $Ne = 0;
	
	foreach my $pos (keys %$pSites) {
		my %bases;
		foreach my $species (@$pSpecies) {
			my $b = substr($pSeq->{$species}, $pos, 1);
			$bases{$b}++;
		}
		my $minorCount = first {defined($_)} sort {$a<=>$b} values %bases;
		my $outB = substr($pSeq->{$outgroup}, $pos, 1);
		
		die if (defined $bases{'N'} or $outB eq 'N');
		
		$N += keys(%bases) - 1;
		$Ns++ if (($minorCount == 1 and not defined $bases{$outB}) or not defined $bases{$outB});
		$Ne++ if ($minorCount == 1 and defined $bases{$outB});
	}
	
	return ($N, $Ns, $Ne);
}

sub getMutationsForSubsampledSites {
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	my $outgroup = shift;
	
	my $N = 0;
	my $Ns = 0;
	my $Ne = 0;
	
	foreach my $pos (keys %$pSites) {
		my $pBases = $pSites->{$pos};
		my $minorCount = first {defined($_)} sort {$a<=>$b} values %$pBases;
		my $outB = substr($pSeq->{$outgroup}, $pos, 1);
		
		die if (defined $pBases->{'N'} or $outB eq 'N');
		
		$N += keys(%$pBases) - 1;
		$Ns++ if (($minorCount == 1 and not defined $pBases->{$outB}) or not defined $pBases->{$outB});
		$Ne++ if ($minorCount == 1 and defined $pBases->{$outB});
	}
	
	return ($N, $Ns, $Ne);
}

sub getPolyWithSubsampling { 
	my $pSeq = shift;
	my $pPoly = shift;
	my $pSpecies = shift;
	my $numSamples = shift;
	
	my $numPoly = 0;

	foreach my $pos (keys %$pPoly) {
		my %bases;
		my $numBases = 0;
		foreach my $s (@$pSpecies) {
			my $b = substr($pSeq->{$s}, $pos, 1);
			if ($b ne 'N') {
				$bases{$b}++;
				$numBases++;
			}
		}
		while ($numBases > $numSamples) {
			my @keys = keys %bases;
			my $b = $keys[int(rand(@keys))];
			$bases{$b}--;
			delete $bases{$b} if ($bases{$b} == 0);
			$numBases--;
		}
		$numPoly++ if (keys %bases > 1);
	}
	
	return $numPoly;
}

sub JukesCantor {
	my $x = shift;
	return (3/4)*log(3/(3-4*$x)); #(-3/4)*log(1-(4*$pi/3))
}

sub Theta {
	my $S = shift;	# number of segregating sites
	my $T = shift;	# total number of sites
	my $n = shift;	# number of sequences
	return unless $T;
	
	my ($harmonic, undef) = harmonic($n);
	return ($S / $harmonic) / $T;
}

sub Pi {  # Pi for all sites
	my $pSeq = shift;
	my $pSpecies = shift;
	
	my $compared = 0;
	my $mismatch = 0;
	
	for (my $i = 0;  $i < @$pSpecies-1;  $i++) {
		my $s1 = $pSpecies->[$i];
		for (my $j = $i+1;  $j < @$pSpecies;  $j++) {
			my $s2 = $pSpecies->[$j];
			die if ($s1 eq $s2);
			for (my $pos = 0;  $pos < length($pSeq->{$s1});  $pos++) {
				my $b1 = substr($pSeq->{$s1}, $pos, 1);
				my $b2 = substr($pSeq->{$s2}, $pos, 1);
				next if ($b1 eq 'N' or $b2 eq 'N'); # skip this site if any missing
				$compared++;
				$mismatch++ if ($b1 ne $b2);
			}
		}
	}

	return 0 if ($compared == 0);
	return $mismatch / $compared;
}

sub PiForSites { # Pi for a list of sites
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	
	my $compared = 0;
	my $mismatch = 0;
	
	for (my $i = 0;  $i < @$pSpecies-1;  $i++) {
		my $s1 = $pSpecies->[$i];
		for (my $j = $i+1;  $j < @$pSpecies;  $j++) {
			my $s2 = $pSpecies->[$j];
			die if ($s1 eq $s2);
			foreach my $pos (keys %$pSites) {
				my $b1 = substr($pSeq->{$s1}, $pos, 1);
				my $b2 = substr($pSeq->{$s2}, $pos, 1);
				next if ($b1 eq 'N' or $b2 eq 'N'); # skip this site if any missing
				$compared++;
				$mismatch++ if ($b1 ne $b2);
			}
		}
	}

	return 0 if ($compared == 0);
	return ($mismatch / $compared);
}

sub PiForSitesWithSubsampling { # Pi for a list of sites with subsampling
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	my $numSamples = shift;
	
	my $compared = 0;
	my $mismatch = 0;

	foreach my $pos (keys %$pSites) {
		my %sample;
		foreach my $s (@$pSpecies) {
			my $b = substr($pSeq->{$s}, $pos, 1);
			$sample{$s} = $b if ($b ne 'N');
		}
		while (keys %sample > $numSamples) {
			my @keys = keys %sample;
			my $i = int(rand(@keys));
			delete $sample{$keys[$i]};
		}
		
		my @species = keys %sample;
		for (my $i = 0;  $i < @species-1;  $i++) {
			my $s1 = $species[$i];
			for (my $j = $i+1;  $j < @species;  $j++) {
				my $s2 = $species[$j];
				die if ($s1 eq $s2);
				$compared++;
				$mismatch++ if ($sample{$s1} ne $sample{$s2});
			}
		}
	}
	
	return 0 if ($compared == 0);
	return ($mismatch / $compared);
}

sub nChoose2 {
	my $n = shift;
	return $n*($n-1)/2;
}

sub PairwiseMismatches { # compute num of pairwise mismatches given bi-allelic polymorphic sites
	my $pPoly = shift;
	my $n = shift;
	
	my $mismatches = 0;
	foreach my $pos (keys %$pPoly) {
		my $minorCount = first {defined($_)} sort {$a<=>$b} values %{$pPoly->{$pos}};
		$mismatches += $minorCount*($n-$minorCount);
	}
	
	return $mismatches;
}

sub PairwiseMismatchesAtSite { # compute num of pairwise mismatches given bi-allelic polymorphic sites
	my $minorCount = shift;
	my $n = shift;
	
	my $mismatches = $minorCount*($n-$minorCount);
	
	return $mismatches;
}

sub PiForPolySites {
	my $pPoly = shift;
	my $len = shift;
	my $n = shift;
	
	my $mismatches = PairwiseMismatches($pPoly, $n);
	my $pairs = nChoose2($n);
	my $comparisons = $len * $pairs;
	return ($mismatches / $comparisons);
}

sub DxyForPolySites { # for bi-allelic polymorphic sites 
	my $pCommonPoly = shift;
	my $pPoly1 = shift;
	my $pPoly2 = shift;
	my $len = shift;
	my $n1 = shift;
	my $n2 = shift;
	
	my $mismatches = 0;
	foreach my $pos (keys %$pCommonPoly) {
		my $b = first {defined($_)} keys %{$pPoly1->{$pos}};
		my $count1 = $pPoly1->{$pos}{$b};
		my $count2 = (defined $pPoly2->{$pos}{$b} ? $pPoly2->{$pos}{$b} : 0);
		my $diff1 = $count1*($n2-$count2);
		my $diff2 = ($n1-$count1)*$count2;
		$mismatches += $diff1 + $diff2;
		#print "pos=$pos $b:$count1,$count2 diff1=$diff1 diff2=$diff2\n";
	}
	
	my $pairs = $n1*$n2;
	my $comparisons = $len * $pairs;
	return ($mismatches / $comparisons);
}

sub DxyForSite { # for bi-allelic polymorphic site
	my $count1 = shift;
	my $count2 = shift;
	my $len = shift;
	my $n1 = shift;
	my $n2 = shift;
	
	my $diff1 = $count1*($n2-$count2);
	my $diff2 = ($n1-$count1)*$count2;
	my $mismatches += $diff1 + $diff2;
	
	my $pairs = $n1*$n2;
	my $comparisons = $len * $pairs;
	return ($mismatches / $comparisons);
}

sub MeanPairwiseDiff { # Mean number of diff sites for a given list of sites
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	
	my $mismatch = 0;
	my $count = 0;
	
	for (my $i = 0;  $i < @$pSpecies-1;  $i++) {
		my $s1 = $pSpecies->[$i];
		for (my $j = $i+1;  $j < @$pSpecies;  $j++) {
			my $s2 = $pSpecies->[$j];
			die if ($s1 eq $s2);
			foreach my $pos (keys %$pSites) {
				my $b1 = substr($pSeq->{$s1}, $pos, 1);
				my $b2 = substr($pSeq->{$s2}, $pos, 1);
				next if ($b1 eq 'N' or $b2 eq 'N'); # skip this site if any missing
				$mismatch++ if ($b1 ne $b2);
			}
			$count++;
		}
	}

	return 0 if ($count == 0);
	return ($mismatch / $count);
}

sub MeanPairwiseDiffWithSubsampling { # for a list of sites with subsampling
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	my $numSamples = shift;
	
	my $count = 0;
	my $mismatch = 0;

	foreach my $pos (keys %$pSites) {
		my %sample;
		foreach my $s (@$pSpecies) {
			my $b = substr($pSeq->{$s}, $pos, 1);
			$sample{$s} = $b if ($b ne 'N');
		} 
		while (keys %sample > $numSamples) {
			my @keys = keys %sample;
			my $i = int(rand(@keys));
			delete $sample{$keys[$i]};
		}
		
		my @species = keys %sample;
		for (my $i = 0;  $i < @species-1;  $i++) {
			my $s1 = $species[$i];
			for (my $j = $i+1;  $j < @species;  $j++) {
				my $s2 = $species[$j];
				die if ($s1 eq $s2);
				$mismatch++ if ($sample{$s1} ne $sample{$s2});
			}
			$count++;
		}
	}
	
	return 0 if ($count == 0);
	return ($mismatch / $count);
}

sub ThetaWithSubsampling { # for a list of sites with subsampling
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies = shift;
	my $numSamples = shift;
	
	my $count = 0;
	my $mismatch = 0;

	foreach my $pos (keys %$pSites) {
		my %sample;
		foreach my $s (@$pSpecies) {
			my $b = substr($pSeq->{$s}, $pos, 1);
			$sample{$s} = $b if ($b ne 'N');
		} 
		while (keys %sample > $numSamples) {
			my @keys = keys %sample;
			my $i = int(rand(@keys));
			delete $sample{$keys[$i]};
		}
		
		my @species = keys %sample;
		for (my $i = 0;  $i < @species-1;  $i++) {
			my $s1 = $species[$i];
			for (my $j = $i+1;  $j < @species;  $j++) {
				my $s2 = $species[$j];
				die if ($s1 eq $s2);
				$mismatch++ if ($sample{$s1} ne $sample{$s2});
			}
			$count++;
		}
	}
	
	return 0 if ($count == 0);
	return ($mismatch / $count);
}

sub DxyForSites { # Average number of pairwise differences per site
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies1 = shift;
	my $pSpecies2 = shift;
	
	my $compared = 0;
	my $mismatch = 0;
	
	foreach my $s1 (@$pSpecies1) {
		foreach my $s2 (@$pSpecies2) {
			die if ($s1 eq $s2);
			foreach my $pos (keys %$pSites) {
				my $b1 = substr($pSeq->{$s1}, $pos, 1);
				my $b2 = substr($pSeq->{$s2}, $pos, 1);
				next if ($b1 eq 'N' or $b2 eq 'N'); # skip this site if any missing
				$compared++;
				$mismatch++ if ($b1 ne $b2);
			}
		}
	}

	return 0 if ($compared == 0);
	return ($mismatch / $compared);
}

sub Fst {
	my $pSeq = shift;
	my $pSpecies1 = shift;
	my $pSpecies2 = shift;
	
	my $within  = Pi($pSeq, $pSpecies1);
	my $total   = Pi($pSeq, [@$pSpecies1, @$pSpecies2]);

	return 0 if ($total == 0);
	my $fst = ($total-$within)/$total;
	return 0 if ($fst <= 0);
	return $fst;
}

sub FstForSites {
	my $pSeq = shift;
	my $pSites = shift;
	my $pSpecies1 = shift;
	my $pSpecies2 = shift;
	
	# Compute weighted average of Pi-within
	my $n = @$pSpecies1 + @$pSpecies2;
	my $w1 = @$pSpecies1 / $n;
	my $w2 = @$pSpecies2 / $n;
	my $within1 = PiForSites($pSeq, $pSites, $pSpecies1);
	my $within2 = PiForSites($pSeq, $pSites, $pSpecies2);
	my $within = $w1*$within1 + $w2*$within2;
	
	# Compute Pi-total
	my $total = PiForSites($pSeq, $pSites, [@$pSpecies1, @$pSpecies2]);

	return 0 if ($total == 0);
	my $fst = ($total-$within)/$total;
	#print "fst=$fst total=$total within=$within within1=$within1 within2=$within2 w1=$w1 w2=$w2\n";
	#return 0 if ($fst <= 0);
	return $fst;
}

sub KaKsForSites { # Ka/Ks for a list of sites
	my $pSites = shift;
	my $pPoly = shift;
	my $pSyn = shift;
	my $pNonsyn = shift;
	
	my ($Pn, $Tn, $Ps, $Ts) = (0, 0, 0, 0);
	
	foreach my $pos (sort {$a<=>$b} keys %$pSites) {
		next if (not defined $pSyn->{$pos});
		
		if (defined $pPoly->{$pos}) {
			if (defined $pNonsyn->{$pos}) { $Pn++; }
			else { $Ps++; }
		}
		$Tn += 1 - $pSyn->{$pos};
		$Ts += $pSyn->{$pos};
	}

	return (0, 0, 0, 0, 0, 0) if ($Tn == 0 or $Ts == 0);
	my $Ka = JukesCantor($Pn/$Tn);
	my $Ks = JukesCantor($Ps/$Ts);
	return return (0, 0, 0, 0, 0, 0) if ($Ks == 0);
	return ($Pn, $Tn, $Ps, $Ts, $Ka, $Ks);
}

sub TajimasD { # Tajima's D 
	my $pi = shift;	# pi value
	my $S = shift;	# number of segregating sites
	my $n = shift;	# number of sequences

	if ($S > 0) {
		# Calculate constants
		my ($a1, $a2) = harmonic($n);
		my $b1 = ($n+1)/(3*($n-1));
		my $b2 = 2*($n*$n + $n + 3)/(9*$n*($n-1));
		my $c1 = $b1 - 1/$a1;
		my $c2 = $b2 - ($n+2)/($a1*$n) + $a2/($a1*$a1);
		my $e1 = $c1/$a1;
		my $e2 = $c2/($a1*$a1 + $a2);
		
		my $theta = $S / $a1;
		my $D = ($pi - $theta)/sqrt($e1*$S + $e2*$S*($S-1));
		return $D;
	}
	
	return 'NaN';
}

sub FuLiD { # Fu & Li's D 
	my $N = shift;
	my $Ns = shift;
	my $Ne = shift;
	my $n = shift; # number of sequences

	if ($N > 0) {
		# Calculate constants
		my ($a, $b) = harmonic($n);
		my ($a2, $b2) = harmonic($n+1);
		my $c = 2*($n*$a - 2*($n - 1)) / (($n-1)*($n-2));
		my $d = $c + (($n-2)/(($n-1)*($n-1))) + (2/($n-1))*(3/2-(2*$a2-3)/($n-2)-1/$n);
		my $v = 1 + ($a*$a/($b+$a*$a))*($c-($n+1)/($n-1));
		my $u = $a - 1 - $v;
		my $v_star = ((($n/($n-1))*($n/($n-1))*$b) + ($a*$a*$d) - (2*($n*$a*($a+1)/(($n-1)*($n-1))))) / ($a*$a + $b);
		my $u_star = ($n/($n-1))*($a - ($n/($n-1))) - $v_star;
			
		my $D = ($N - $a*$Ne) / sqrt($u*$N + $v*$N*$N);
		my $D_star = (($n/($n-1))*$N - $a*$Ns) / sqrt($u_star*$N + $v_star*$N*$N);
		my $sigD = testSignificance($D, $n, \%D_significance);
		my $sigDStar = testSignificance($D_star, $n, \%DStar_significance);
		return ($D, $sigD, $D_star, $sigDStar);
	}
	
	return ('NaN', 0, 'NaN', 0);
}

sub testSignificance {
	my $val = shift;
	my $n = shift;
	my $pTable = shift;
	
	die if (not defined $pTable->{$n});
	
	my $upper = $pTable->{$n}->[3];
	my $lower = $pTable->{$n}->[2];
	return ($val < $lower or $val > $upper);
}

sub harmonic {
	my $n = shift;

	my ($h, $h2) = (0, 0);
	for (my $i = 1;  $i < $n;  $i++) {
		$h += 1/$i;
		$h2 += 1/($i*$i);
	}
	
	return ($h, $h2);
}

1;