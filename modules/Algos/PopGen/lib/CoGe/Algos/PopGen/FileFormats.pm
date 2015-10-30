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

our @EXPORT = qw( loadGFF loadBED loadCounts loadVCF parseVCF );

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

		$chr = _formatChr($chr);
		
		# Parse attributes field
		my %attr = $attr =~ /(\w+)=([^;]+);*/g;
		%attr = map { lc($_) => $attr{$_} } keys %attr;
		my $id = $attr{id};
		my $parent = $attr{parent};
		$parent = $id unless $parent;
		
		#print "Warning: duplicate gene id $id for '$name'\n" if (defined $out{$type}{$id}) if ($type eq 'gene');
		push @{$out{$type}{$parent}}, { 
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

sub loadVCF {
    my $filename = shift;
    
    open(my $fh, $filename) or 
        die("Error: cannot open file '$filename'\n");
    
    my (%variants, %sites);
    while (my $line = <$fh>) {
        next if ($line =~ /^#/);
        chomp $line;
        my ($chr, $pos, undef, $ref, $alt_str, undef, undef, undef, $format_str, $genotype_str) = split(/\t/, $line, 10);
        $chr = _formatChr($chr);
        
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
    $chr = _formatChr($chr);

    my @genotypes;
    my @format = split(':', $format_str);
    foreach (split(/\t/, $genotype_str)) {
        my $i = 0;
        my %gt = map { $format[$i++] => $_ } split(':');
        push @genotypes, \%gt;
    }
    #print STDERR Dumper \@genotypes, "\n";
    
    my %counts;
    my @alleles = split(',', $alt_str);
    unshift @alleles, $ref;
    foreach my $gt (@genotypes) {
        next unless ($gt->{GT} && $gt->{GT} ne './.');
        my @altIndex = split('/', $gt->{GT});
        map { $counts{ $alleles[$_] }++  } @altIndex;
    }
    #print STDERR $line, "\n", Dumper \%counts, "\n";

    return {
        'chr' => $chr,
        'pos' => $pos,
        'ref' => $ref,
        'alleles' => \%counts
    };
}

sub _formatChr {
    my $chr = shift;
    $chr = lc($chr);
    #$chr =~ s/_//;
    $chr =~ s/(\d+)/\_$1/;
    return $chr;   
}

1;