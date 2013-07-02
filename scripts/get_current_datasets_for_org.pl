#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $coge = CoGeX->dbconnect();

my $oid = shift;

my ($org) = $coge->resultset('Organism')->resolve($oid);
print $org->name;
print ": ",$org->description if $org->description;
print "\n";
my @items;
foreach my $ds ($org->current_datasets())
  {
    my $str;
    $str .= join (", ", $ds->get_chromosomes);

    $str .=  "\tv:".$ds->version. "\t". $ds->name;
    $str.=  ": ".$ds->description if $ds->description;
    $str .= "\n";
    push @items, $str if $str;
  }
print sort @items,"\n";
