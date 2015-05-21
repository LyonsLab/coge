#!/usr/bin/perl -w

use strict;
use CoGeX;

my $org = shift;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

($org) = $coge->resultset('Organism')->resolve($org);

foreach my $ds ($org->current_datasets)
  {
    my $name = $ds->name;
    $name .= ": ".$ds->description if $ds->description;
    $name = "<a href=".$ds->link." target=_new>".$name."</a>"if $ds->link;
    my $source = $ds->datasource->name;
    $source .= ": ".$ds->datasource->description if $ds->datasource->description;
    $source = "<a href=".$ds->datasource->link."target=_new>".$source."</a>" if $ds->datasource->link;
    print "Current datasets for ".$org->name;
    print ": ".$org->description if $org->description;
    print "<br>\n";
    print join ("\t", $name, $source),"\n";
  }
