#!/usr/bin/perl -w
#Eric Lyons 2017
#Syntenic Path Assembly algorithm
#Uses synteny output from DagChainer to generate a syntenic path assembly

use strict;
use Getopt::Long;
use CoGe::Accessory::SynMap_report;
use Data::Dumper;

use vars qw($outputfile $alignfile $reference $skip_random);

GetOptions(
	"output=s"			=>	\$outputfile,
	"alignfile=s"			=>	\$alignfile,
        #"skip_non_aligned_contigs"	=>	\$skip_non_aligned_contigs,
        "reference"                     =>      \$reference, #default 0 will align genome with more pieces to genome with fewer pieces.  "1" aligns the other way
,
	"skip_random"			=>	\$skip_random,
);

#defaults
#$skip_non_aligned_contigs = 0 unless $skip_non_aligned_contigs;
$reference = 0 unless $reference;
$skip_random = 0 unless $skip_random;

unless (-r $alignfile) {
  usage();
}

my $synmap_report = new CoGe::Accessory::SynMap_report;

    print STDERR 'alignfile=', $alignfile, "\n";
    my ( $org1_association, $org2_association );
    ( $org1_association, $org2_association ) =
        $synmap_report->parse_syn_blocks( file => $alignfile );
    my $output;
    if (   (@$org1_association > @$org2_association && $reference == 0)
        || (@$org1_association < @$org2_association && $reference == 1) )
    {
        $output = reord(
            assoc       => $org2_association,
            info        => $org1_association,
            skip_random => $skip_random
        );
    }
    else {
        $output = reord(
            assoc       => $org1_association,
            info        => $org2_association,
            skip_random => $skip_random
        );
    }

    if ($outputfile) {
      open OUT, ">$outputfile";
      print OUT $output;
      close OUT;
    }
    else {
      print $output;
    }


sub reord {
    my %opts        = @_;
    my $association = $opts{assoc};
    #my $skip        = $opts{skip};
    my $skip_random = $opts{skip_random};
    my $info        = $opts{info};          #for determining orientation
    #print STDERR 'reord: info=', scalar(@$info), ' assoc=', scalar(@$association), ' skip=', $skip, ' skip_random=', $skip_random, "\n";

    #create mapping hash of which contigs are in the reverse orientation
    my %rev_info;
    foreach my $item (@$info) {
        my $rev = $item->{rev} ? -1 : 1;
        my $chr = $item->{chr};
        $rev_info{$chr} = $rev;
    }

    my %mapped_association;
    map { $mapped_association{ $_->{chr} } = $_->{matching_chr} } @$association;
    my @new_order;

    my $output = join( "\t", ( "#CHR1", "CHR2", "ORIENTATION" ) ) . "\n";
    foreach my $chr (map ($_->{chr}, @$association)) {
        next if check_random($chr) && $skip_random;
        if ( $mapped_association{$chr} ) {
            foreach my $item ( @{ $mapped_association{$chr} } ) {
                next if check_random($item) && $skip_random;
            }
        }
        $output .= join( "\n",
            map { join( "\t", $chr, $_, $rev_info{$_} ) }
              @{ $mapped_association{$chr} } )
          . "\n";
    }
    #unless ($skip) {
    #    my %seen = map { $_ => 1 } @new_order;
    #    foreach my $item (@$reorder) {
    #        $output .= join( "\t", "unmapped", $item ) . "\n"
    #          unless $seen{$item};
    #        push @new_order, $item unless $seen{$item};
    #    }
    #}
    #@$reorder = @new_order;
    return $output;
}

sub check_random {
    my $chr = shift;
    return 1 if ( $chr =~ /random/i || $chr =~ /unknown/i || $chr =~ /^un$/i );
    return 0;
}

sub usage {
 print qq{Usage: $0 -alignfile <synmap alignment file>

Description:
 Syntenic Path Assembly algorithm
 Uses synteny output from DagChainer to generate a syntenic path assembly
 Details on algorithm and implementation in CoGe: https://genomevolution.org/wiki/index.php/Syntenic_path_assembly

Options:
 -alignfile <syntenic alignment file           :  syntenic file from Synmap, usually named: XXXX.aligncoords.gcoords
 -output <output file>                         :  send output to a file.  If not used, will print to STDOUT
 -reference                                    :  which genome will be used as the reference genome.  
						  0 / undef:  (default) will use the genome with fewer pieces/chromosomes
                                                  1:          Use the genome with more pieces/chromosomes
 -skip_random				       :  Skip chromosomes/contigs that contain the work "random" or "unknown"
                                                  0 / undef: (default) don't skip those chromosomes
						  1:         Skip those chromosomes (will not be included in output)
Limitations:
 Chromosomes that don't have syntenic matches are not listed
 Order of chromosomes from reference genome is random (no sorting).
};
 exit(0);
}


