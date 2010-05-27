#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

foreach my $feat ($coge->resultset('FeatureType')->find(185)->features())
  {
    $feat->delete;
  }
