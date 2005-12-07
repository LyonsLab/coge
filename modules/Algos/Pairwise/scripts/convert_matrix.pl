#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my %data;
my @head;
while (<>)
  {
    next if /^#/;
    chomp;
    next unless $_;
    my @line = split /\s+/;
    my $chr = shift @line;
    unless ($chr)
      {
	@head = @line;
	next;
      }
    my $i=0;
    foreach my $score (@line)
      {
	$data{$chr}{$head[$i]} = $score;
	$i++;
      }
  }
print Dumper \%data;
