#!/usr/bin/perl -w

use strict;
use CoGeX;

my $dsid = shift;
unless ($dsid)
  {
    print qq{
Usage: $0 <dataset name or id>

This program will remove the specified dataset and all related data from the CoGe genomes database.
};
    exit;
  }

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
my ($ds) = $coge->resultset('Dataset')->resolve($dsid);
unless ($ds)
  {
    print "Unable to find a dataset for $dsid\n";
    exit;
  }
print "Do you want to delete:  \n\t".$ds->name.": ".$ds->description." from organism".$ds->organism->name,": ".$ds->organism->description," (y/n)?\n";
my $var = <STDIN>;
if ($var =~ /y/i)
  {
    print "Deleteing. . .";
    $ds->delete;
    print "completed!\n";
  }
else
  {
    print "Delete aborted.\n";
  }
