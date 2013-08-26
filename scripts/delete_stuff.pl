#!/usr/bin/perl -w

use strict;

use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $org = $coge->resultset('Organism')->find('4310');
foreach my $ds ($org->datasets)
  {
    $ds->delete if $ds->name =~ /^NW_/;
  }
