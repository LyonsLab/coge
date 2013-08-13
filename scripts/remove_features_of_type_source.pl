#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

foreach my $feat ($coge->resultset('FeatureType')->find(185)->features())
  {
    $feat->delete;
  }
