#-------------------------------------------------------------------------------
# Purpose:	Functions for loading various file types
# Author:	Matt Bomhoff
# Created:	10/9/15
#-------------------------------------------------------------------------------
package CoGe::Algos::PopGen::FileFormats;

use warnings;
use strict;
use Data::Dumper;
#use Set::IntSpan::Fast;
use base 'Exporter';

our @EXPORT = qw( loadFASTA loadGFF loadBED loadCounts detectGVCF loadVCF parseVCF parseGenotypes );

sub loadFASTA {
    my $filename = shift;
    my %out;
        
    open(my $fh, $filename) or 
        die("Error: cannot open file '$filename'\n");
    
    my $name;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^\>(\S+)/) { # FASTA header
            $name = $1;
            die "Error: sequences with same name in input file\n" if (defined $out{$name});
        }
        else { # FASTA data
            die "loadFASTA: parse error" if (not defined $name);
            my ($s) = $line =~ /\S+/g; # because of Windows end-of-line
            $out{$name} .= $s if (defined $s);
        }
    }
    close($fh);
    
    return \%out;
}

sub loadGFF {
	my %opts = @_;
	my $file       = $opts{file};
	
	open(my $fh, $file) or 
		die("Error: cannot open file '$file'\n");
	
	my %out;
	while (<$fh>) {
	    next if /^#/;
		chomp;
		my ($chr, undef, $type, $start, $end, undef, $strand, undef, $attr) = split /\t/;
		next unless (defined $chr && defined $type && defined $start && defined $end && defined $attr);
		$type = lc($type);
		($start, $end) = ($end, $start) if ($end < $start); # fix coordinates if needed

#		$chr = _formatChr($chr);
		
		# Parse attributes field
		my %attr = $attr =~ /(\w+)=([^;]+);*/g;
		%attr = map { lc($_) => $attr{$_} } keys %attr;
		my $id = $attr{id};
		my $parent = $attr{parent};
		$parent = $id unless $parent;
		
		#print "Warning: duplicate gene id $id for '$name'\n" if (defined $out{$type}{$id}) if ($type eq 'gene');
		push @{$out{$type}{$chr}{$parent}}, { 
		    'chr'    => $chr, 
		    'start'  => $start, 
		    'end'    => $end, 
		    'strand' => $strand, 
		    'id'     => $id, 
		    'attr'   => \%attr
		};
	}
	close($fh);
	
#	print "" . (keys %out) . " total genes loaded\n";

	# Remove where CDS w/o gene
#	foreach my $name (keys %out) {
#		if (not defined $out{$name}{'chr'}) {
#			delete $out{$name};
#		}
#	}
	
	return \%out;
}

sub loadBED {
    my $filename = shift;
    
    open(my $fh, $filename) or 
        die("Error: cannot open file '$filename'\n");
    
    my %out;
    while (<$fh>) {
        next if /^#/;
        chomp;
        my ($chr, $start, $end) = split(/\t/);
        $chr = _formatChr($chr);
        
        foreach ($start+1..$end) {
            $out{$chr}{$_} = 1;
        }
    }
    
    close($fh);
    
    return \%out;
}

sub loadCounts {
    my $filename = shift;
    
    open(my $fh, $filename) or 
        die("Error: cannot open file '$filename'\n");
    
    my %out;
    <$fh>; # skip header line
    while (<$fh>) {
        chomp;
        my ($chr, $pos, undef, undef, $allele_str) = split(/\t/, $_, 5);
        $chr = _formatChr($chr);
        
        my %alleles;
        foreach (split(/\t/, $allele_str)) {
            my ($base, $count) = split(':');
            $alleles{$base} = $count;
        }
        $out{$chr}{$pos} = \%alleles;
    }
    
    close($fh);
    
    return \%out;
}

sub detectGVCF {
    my $filename = shift;
    my $NUM_LINES_TO_SAMPLE = 100;

    open(my $fh, $filename) or
        die("Error: cannot open file '$filename'\n");

    my ($count, $isGVCF);
    while (my $line = <$fh>) {
        next if ($line =~ /^#/);
        chomp $line;
        my @tok = split(/\t/, $line);
        $isGVCF = 1 if (@tok > 10); # this is a GVCF file, more than one sample detected
        goto DONE if ($isGVCF || $count++ > $NUM_LINES_TO_SAMPLE);
    }
    DONE:
    close($fh);
    return $isGVCF;
}

sub loadVCF {
    my $filename = shift;
    
    open(my $fh, $filename) or 
        die("Error: cannot open file '$filename'\n");
    
    my (%variants, %sites);
    while (my $line = <$fh>) {
        next if ($line =~ /^#/);
        chomp $line;
        my ($chr, $pos, undef, $ref, $alt_str, undef, undef, undef, $format_str, $genotype_str) = split(/\t/, $line, 10);
#        $chr = _formatChr($chr);
        
        $sites{$chr}{$pos}++;
        next if ($alt_str eq '.');

        $variants{$chr}{$pos} = {
            'ref'    => $ref,
            'alt'    => $alt_str,
            'format' => $format_str,
            'gt'     => $genotype_str
        };
    }
    
    close($fh);
    
    return (\%variants, \%sites);
}

sub parseVCF {
    my $line = shift; # VCF string (one line)
    
    chomp $line;
    my ($chr, $pos, undef, $ref, $alt_str, undef, undef, undef, $format_str, $genotype_str) = split(/\t/, $line, 10);
    #$chr = _formatChr($chr);

    return {
        'chr' => $chr,
        'pos' => $pos,
        'ref' => $ref,
        'alt' => $alt_str,
        'format' => $format_str,
        'genotypes' => $genotype_str
    };
}

sub parseGenotypes {
    my $vcf = shift; # parsed VCF record
    
    my @genotypes;
    my @format = split(':', $vcf->{format});
    foreach (split(/\t/, $vcf->{genotypes})) {
        my $i = 0;
        my %gt = map { $format[$i++] => $_ } split(':');
        push @genotypes, \%gt;
    }
    #print STDERR Dumper \@genotypes, "\n";
    
    my %counts;
    my @alleles = split(',', $vcf->{alt});
    unshift @alleles, $vcf->{'ref'};
    foreach my $gt (@genotypes) {
        next unless ($gt->{GT} && $gt->{GT} ne './.');
        my @altIndex = split('/', $gt->{GT});
        map { $counts{ $alleles[$_] }++  } @altIndex;
    }
    #print STDERR $line, "\n", Dumper \%counts, "\n";
    
    return \%counts;
}

sub _formatChr {
    my $chr = shift;
    $chr = lc($chr);
    #$chr =~ s/_//;
    $chr =~ s/(\d+)/\_$1/; # for naming inconsistency between Cgrandiflora_266_v1.1.gene and Cg_scattered.vcf
    return $chr;   
}

1;