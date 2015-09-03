#!/usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';

use Getopt::Long;
use CoGe::Accessory::Web;

our ($infile, $outfile, $max_distance, $P, $QUOTA_ALIGN, $CLUSTER_UTILS,
    $CONFIG);

GetOptions(
    "infile|i=s"       => \$infile,
    "outfile|o=s"      => \$outfile,
    "max_distance|d=s" => \$max_distance,
    "config|cfg=s"     => \$CONFIG,);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
$ENV{PATH} = join ":",
  (
    $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin");

$QUOTA_ALIGN   = $P->{QUOTA_ALIGN};     #the program
$CLUSTER_UTILS = $P->{CLUSTER_UTILS};   #convert dag output to quota_align input

run_quota_align_merge(
    infile   => $infile,
    outfile  => $outfile,
    max_dist => $max_distance);

sub run_quota_align_merge
{
    my %opts     = @_;
    my $infile   = $opts{infile};
    my $max_dist = $opts{max_dist};
    my $outfile  = $opts{outfile};

    #convert to quota-align format
    my $cmd = $CLUSTER_UTILS
      . " --format=dag --log_evalue $infile $infile.Dm$max_dist.qa";
    `$cmd`;

    say "Converting dag output to quota_align format: $cmd";
    $cmd = $QUOTA_ALIGN . " --Dm=$max_dist --merge $infile.Dm$max_dist.qa";

    say "Running quota_align to merge diagonals:\n\t$cmd";
    `$cmd`;

    if (-r "$infile.Dm$max_dist.qa.merged") {
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
        open(IN,  "$infile.Dm$max_dist.qa.merged");
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
    } else {
        say "The merged file $infile.Dm$max_dist.qa.merged was not created.";
    }
}
