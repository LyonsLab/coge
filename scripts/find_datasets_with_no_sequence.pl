#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my $coge = CoGeX->dbconnect();

foreach my $ds ($coge->resultset('Dataset')->all())
  {
    my @chrs = $ds->get_chromosomes();
    unless (@chrs)
      {
	print $ds->organism->name," (",$ds->organism->id,")",": ",$ds->name," (",$ds->id,")\n";
      }
    foreach my $chr (@chrs)
      {
	my $last = $ds->last_chromosome_position($chr);
	next if $last;
	print $ds->organism->name,": ",$ds->name,"\n";
      }
  }
