#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $rs = $coge->resultset('Feature');
while( my $feat=$rs->next)
  {
    foreach my $loc ($feat->locations)
      {
	print "strand mismatch ",$feat->strand,"::",$loc->strand,"\n" if $feat->strand ne $loc->strand;
      }
  }
