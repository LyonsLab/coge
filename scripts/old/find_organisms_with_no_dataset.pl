#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

foreach my $org ($coge->resultset('Organism')->all())
  {
    my @ds =  $org->datasets;
    next if scalar @ds;
    print "Deleting: ", $org->name,": ",$org->description,"\n";
#    $org->delete();
  }
