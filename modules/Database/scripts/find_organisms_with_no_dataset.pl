#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

foreach my $org ($coge->resultset('Organism')->all())
  {
    my @ds =  $org->datasets;
    next if scalar @ds;
    print "Deleting: ", $org->name,": ",$org->description,"\n";
#    $org->delete();
  }
