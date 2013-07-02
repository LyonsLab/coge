#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=coge;host=biocon.berkeley.edu;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

foreach my $org ($coge->resultset('Organism')->all())
  {
    my @dsg =  $org->dataset_groups;
    next if scalar @dsg;
    print "Deleting: ", $org->name,": ",$org->description,"\n";
    $org->delete();
  }
