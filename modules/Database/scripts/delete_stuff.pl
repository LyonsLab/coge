#!/usr/bin/perl -w

use strict;

use CoGeX;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $org = $coge->resultset('Organism')->find('4310');
foreach my $ds ($org->datasets)
  {
    $ds->delete if $ds->name =~ /^NW_/;
  }
