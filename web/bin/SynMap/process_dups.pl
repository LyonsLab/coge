#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

our ($cogeweb, $input, $output, $TEMPDIR, $CONFIG, $P, $BASE_URL);

use Getopt::Long;
use CoGe::Accessory::Web;

GetOptions(
    "infile=s"     => \$input,
    "outfile=s"     => \$output,
    "config|cfg=s" => \$CONFIG,
);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
$TEMPDIR = $P->{TEMPDIR} . "SynMap";
$BASE_URL = $P->{SERVER};

$cogeweb = CoGe::Accessory::Web::initialize_basefile(
    tempdir  => $TEMPDIR
);
$| = 1;    # Disable buffering

process_local_dups_file(infile => $input, outfile => $output);

sub process_local_dups_file {
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};

    CoGe::Accessory::Web::write_log(
"Adding coge links to tandem duplication file.  Infile $infile : Outfile $outfile",
        $cogeweb->logfile
    );
    $/ = "\n";
    open( IN,  $infile );
    open( OUT, ">$outfile" );
    print OUT "#",
      join( "\t",
        "FeatList_link", "GEvo_link", "FastaView_link",
        "chr||start||stop||name||strand||type||database_id||gene_order" ),
      "\n";

    while (<IN>) {
        chomp;
        next unless $_;
        my @line = split /\t/;
        my %fids;
        foreach (@line) {
            my @item = split /\|\|/;
            next unless $item[6];
            $fids{ $item[6] } = 1;
        }
        next unless keys %fids;
        my $featlist = $BASE_URL . "FeatList.pl?";
        map { $featlist .= "fid=$_;" } keys %fids;
        my $fastaview = $BASE_URL . "FastaView.pl?";
        map { $fastaview .= "fid=$_;" } keys %fids;
        my $gevo  = $BASE_URL . "GEvo.pl?";
        my $count = 1;

        foreach my $id ( keys %fids ) {
            $gevo .= "fid$count=$id;";
            $count++;
        }
        $gevo .= "num_seqs=" . scalar keys %fids;
        print OUT join( "\t", $featlist, $gevo, $fastaview, @line ), "\n";
    }
    close OUT;
    close IN;
}
