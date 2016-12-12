#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

use CoGe::Accessory::Web;
use Getopt::Long;

our ( $cogeweb, $input, $output, $TEMPDIR, $CONFIG, $P, $BASE_URL, $id1, $id2 );

GetOptions(
    "infile|i=s"   => \$input,
    "outfile|o=s"  => \$output,
    "dsgid1|id1=s" => \$id1,
    "dsgid2|id2=s" => \$id2,
    "config|cfg=s" => \$CONFIG,
);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
unless ($P) {
    print STDERR "Couldn't load config file\n";
    exit(-1);
}

unless ($id1 && $id2) {
    print STDERR "Missing genome IDs\n";
    exit(-1);    
}

unless ($input && $output && -e $input) {
    print STDERR "Missing input/output path\n";
    exit(-1);      
}

$TEMPDIR  = $P->{TEMPDIR} . "SynMap";
$BASE_URL = $P->{SERVER};

$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
$| = 1;    # Disable buffering

generate_links(
    infile  => $input,
    outfile => $output,
    dsgid1  => $id1,
    dsgid2  => $id2
);

sub generate_links {
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};
    my $dsgid1  = $opts{dsgid1};
    my $dsgid2  = $opts{dsgid2};
    $/ = "\n";

    #
    # Condense hits and add GEvo links to input file
    #

    open( IN,  $infile );
    open( OUT, ">$outfile" );
    my %condensed;
    my %names;

    while (<IN>) {
        chomp;
        if (/^#/) {
            print OUT $_, "\n";
            next;
        }

        s/^\s+//;
        next unless $_;
        my @line  = split /\t/;
        my @feat1 = split /\|\|/, $line[1];
        my @feat2 = split /\|\|/, $line[5];
        my $link  = $BASE_URL . "GEvo.pl?";
        my ( $fid1, $fid2 );

        if ( $feat1[6] ) {
            $fid1 = $feat1[6];
            $link .= "fid1=" . $fid1;
        }
        else {
            my ($xmin) = sort ( $feat1[1], $feat1[2] );
            my $x = sprintf( "%.0f", $xmin + abs( $feat1[1] - $feat1[2] ) / 2 );
            $link .= "chr1=" . $feat1[0] . ";x1=" . $x;
        }
        if ( $feat2[6] ) {
            $fid2 = $feat2[6];
            $link .= ";fid2=" . $fid2;
        }
        else {
            my ($xmin) = sort ( $feat2[1], $feat2[2] );
            my $x = sprintf( "%.0f", $xmin + abs( $feat2[1] - $feat2[2] ) / 2 );
            $link .= ";chr2=" . $feat2[0] . ";x2=" . $x;
        }
        $link .= ";dsgid1=" . $dsgid1;
        $link .= ";dsgid2=" . $dsgid2;

        if ( $fid1 && $fid2 ) {
            $condensed{ $fid1 . "_" . $dsgid1 }{ $fid2 . "_" . $dsgid2 } = 1;
            $condensed{ $fid2 . "_" . $dsgid2 }{ $fid1 . "_" . $dsgid1 } = 1;
            $names{$fid1} = $feat1[3];
            $names{$fid2} = $feat2[3];
        }

        #accn1=".$feat1[3]."&fid1=".$feat1[6]."&accn2=".$feat2[3]."&fid2=".$feat2[6] if $feat1[3] && $feat1[6] && $feat2[3] && $feat2[6];
        print OUT $_;
        print OUT "\t", $link;
        print OUT "\n";
    }
    close IN;
    close OUT;

    #
    # Generate file of condensed hits
    #

    if ( #keys %condensed &&
        !( -r "$outfile.condensed" || -r "$outfile.condensed.gz" ) )
    {
        open( OUT, ">$outfile.condensed" );
        print OUT
            join( "\t",
            qw(COUNT GEVO MASKED_GEVO FASTA_LINK GENE_LIST GENE_NAMES) ),
            "\n";

        #take into account transitivity
        foreach my $id2 ( keys %condensed ) {
            foreach my $id2 ( keys %{ $condensed{$id1} } ) {
                foreach my $id3 ( keys %{ $condensed{$id2} } ) {
                    next if $id1 eq $id2;
                    $condensed{$id1}{$id3} = 1;
                    $condensed{$id3}{$id1} = 1;
                }
            }
        }

        my %seen;
        foreach my $id1 (
            sort {
                scalar( keys %{ $condensed{$b} } ) <=> scalar( keys %{ $condensed{$a} } )
            } keys %condensed
          )
        {
            my ( $fid1, $dsgid1 ) = split /_/, $id1;
            next if $seen{$fid1};
            $seen{$fid1} = 1;
            my @names     = $names{$fid1};
            my $gevo_link = $BASE_URL . "GEvo.pl?fid1=$fid1;dsgid1=$dsgid1";
            my $fids      = "fid=$fid1";
            my $count     = 2;

            foreach my $id2 ( sort keys %{ $condensed{$id1} } ) {
                my ( $fid2, $dsgid2 ) = split /_/, $id2, 2;
                next if $fid1 == $fid2;
                $seen{$fid2} = 1;
                $gevo_link .= ";fid$count=$fid2;dsgid$count=$dsgid2";
                $fids      .= ",$fid2";
                push @names, $names{$fid2};
                $count++;
            }
            $count--;
            $gevo_link .= ";num_seqs=$count";
            my $gevo_link2 = $gevo_link;
            $gevo_link .= ";pad_gs=20000";
            for my $i ( 1 .. $count ) {
                $gevo_link2 .= ";mask$i=non-cds";
            }
            $gevo_link2 .= ";pad_gs=200000";
            $gevo_link2 .= ";autogo=1";
            my $fasta_link    = $BASE_URL . "FastaView.pl?$fids";
            my $featlist_link = $BASE_URL . "FeatList.pl?$fids";
            print OUT join( "\t",
                $count, $gevo_link, $gevo_link2, $fasta_link, $featlist_link,
                @names ),
                "\n";
        }
    }
    close OUT;
}
