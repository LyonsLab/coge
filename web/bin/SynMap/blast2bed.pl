#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $infile;
my $outfile1;
my $outfile2;

GetOptions ( 
	    "infile=s" => \$infile,
	    "outfile1=s"    => \$outfile1,
	    "outfile2=s"    => \$outfile2,
	   );

open (OUT1, ">$outfile1") || die "Can't open $outfile1 for writing: $!";
open (OUT2, ">$outfile2") || die "Can't open $outfile2 for writing: $!";
open (IN, $infile)  || die "Can't open $infile for reading: $!";
my %seen1 = ();
my %seen2 = ();
while (<IN>)
  {
    chomp;
    next unless $_;
    my @line = split/\t/;
    my @item1 = split/\|\|/, $line[0];
    my @item2 = split/\|\|/, $line[1];
    #genomic comparisons won't have starts and stops in the name file, replace those with the actual hits
    $item1[1] = $line[6] unless defined $item1[1];
    $item1[2] = $line[7] unless defined $item1[2];
    $item2[1] = $line[8] unless defined $item2[1];
    $item2[2] = $line[9] unless defined $item2[2];
    print OUT1 join ("\t", $item1[0], $item1[1], $item1[2], $line[0]),"\n" unless $seen1{$line[0]};
    print OUT2 join ("\t", $item2[0], $item2[1], $item2[2], $line[1]),"\n" unless $seen2{$line[1]};
    $seen1{$line[0]} =1;
    $seen2{$line[1]} =1;
  }
close IN;
close OUT1;
close OUT2;
