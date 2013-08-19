#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $rs = $coge->resultset('Feature');
while( my $feat=$rs->next)
  {
    foreach my $loc ($feat->locations)
      {
	print "strand mismatch ",$feat->strand,"::",$loc->strand,"\n" if $feat->strand ne $loc->strand;
      }
  }
