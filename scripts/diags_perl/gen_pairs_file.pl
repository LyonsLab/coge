#!/usr/bin/perl -w

use strict;
use CoGeX;

my $coge = CoGeX->dbconnect();

print join ("\t", qw(CHR1 START1 STOP1 NAME1 CHR2 START2 STOP2 NAME2)),"\n";

while (<>)
  {
    next if /^#/;
    chomp;
    next unless $_;
    my @line = split /\t/;
    my @item1 = split /\|\|/,$line[1];
    my @item2 = split /\|\|/,$line[5];
    my $feat1 = $coge->resultset('Feature')->find($item1[6]);
    my $feat2 = $coge->resultset('Feature')->find($item2[6]);
    my $name1;
    foreach ($feat1->names)
      {
	$name1 = $_ unless $name1;
	next if /_T/;
	$name1 = $_ if length ($_) > length $name1;
      }
    my $name2;
    foreach ($feat2->names)
      {
	$name2 = $_ unless $name2;
	next if /_T/;
	$name2 = $_ if length ($_) > length $name2;
      }
    print join ("\t", $feat1->chromosome, $feat1->start, $feat1->stop, $name1, $feat2->chromosome, $feat2->start, $feat2->stop, $name2),"\n";
  }
