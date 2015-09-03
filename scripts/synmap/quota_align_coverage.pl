#!/usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';

use Getopt::Long;
use CoGe::Accessory::Web;

our (
    $infile,      $outfile,       $orgratio1,
    $orgratio2,   $overlap,       $P,
    $QUOTA_ALIGN, $CLUSTER_UTILS, $CONFIG);

GetOptions(
    "infile|if=s"           => \$infile,
    "outfile|of=s"          => \$outfile,
    "depth_ratio_org1|d1=s" => \$orgratio1,
    "depth_ratio_org2|d2=s" => \$orgratio2,
    "depth_overlap|o=s"     => \$overlap,
    "config|cfg=s"          => \$CONFIG,);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
$ENV{PATH} = join ":",
  (
    $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin");

$QUOTA_ALIGN   = $P->{QUOTA_ALIGN};     #the program
$CLUSTER_UTILS = $P->{CLUSTER_UTILS};   #convert dag output to quota_align input

run_quota_align_coverage(
    infile  => $infile,
    outfile => $outfile,
    org1    => $orgratio1,
    org2    => $orgratio2,
    overlap => $overlap);

sub run_quota_align_coverage
{
    my %opts         = @_;
    my $infile       = $opts{infile};
    my $org1         = $opts{org1};      #ratio of org1
    my $org2         = $opts{org2};      #ratio of org2
    my $overlap_dist = $opts{overlap};
    my $outfile      = $opts{outfile};

    #convert to quota-align format
    my $cov_cmd =
      $CLUSTER_UTILS . " --format=dag --log_evalue $infile $infile.qa";
    my $qa_cmd = $QUOTA_ALIGN
      . " --Nm=$overlap_dist --quota=$org1:$org2 $infile.qa > $outfile.tmp";

    say "Convert command: $cov_cmd";
    say "Quota Align command: $qa_cmd";

    say "Converting dag output to quota_align format.";
    `$cov_cmd`;

    say "Running quota_align to find syntenic coverage.";
    my $qa_output = `$qa_cmd`;

    if (-r "$outfile.tmp") {
        my %data;
        $/ = "\n";
        open(IN, $infile);
        while (<IN>) {
            next if /^#/;
            my @line = split /\t/;
            $data{join("_", $line[0], $line[2], $line[4], $line[6])} = $_;
        }
        close IN;
        open(OUT, ">$outfile");
        open(IN,  "$outfile.tmp");
        while (<IN>) {
            if (/^#/) {
                print OUT $_;
            } else {
                chomp;
                my @line = split /\t/;
                print OUT $data{
                    join("_", $line[0], $line[1], $line[2], $line[3])};
            }
        }
        close IN;
        close OUT;

        generate_grimm_input(infile => $outfile);
    } else {
        say "Syntenic coverage failed to output $outfile.tmp";
    }
}

sub generate_grimm_input
{
    my %opts   = @_;
    my $infile = $opts{infile};

    CoGe::Accessory::Web::gunzip($infile . ".gz")    if -r $infile . ".gz";
    CoGe::Accessory::Web::gunzip($infile . ".qa.gz") if -r $infile . ".qa.gz";
    my $cmd = $CLUSTER_UTILS . " --format=dag --log_evalue $infile $infile.qa";
    say "\nGenerating input data for GRIMM";
    say "Converting dag output to quota_align format: $cmd";
    `$cmd`;
    $cmd = $CLUSTER_UTILS . " --print_grimm $infile.qa";
    say "running  cluster_utils to generating grimm input:\n\t$cmd";
    my $output;
    open(IN, "$cmd |");

    while (<IN>) {
        $output .= $_;
    }
    close IN;
    my @seqs;
    foreach my $item (split /\n>/, $output) {
        $item =~ s/>//g;
        my ($name, $seq) = split /\n/, $item, 2;
        $seq =~ s/\n$//;
        push @seqs, $seq;
    }

    return \@seqs;
}
