#!/usr/bin/perl -w

use strict;
use CoGeX;

my $orgid = shift;
unless ($orgid)
  {
    print qq{
Usage: $0 <organism name or id>

This program will remove the specified organism and all related data from the CoGe genomes database.
};
    exit;
  }

my $connstr = 'dbi:mysql:dbname=coge;host=homer;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
my $org = $coge->resultset('Organism')->resolve($orgid);

print "Do you want to delete:  \n\t".$org->name.": ".$org->description."?\n";
my $var = <STDIN>;
if ($var =~ /y/i)
  {
    print "Deleteing. . .";
    $org->delete;
    print "completed!\n";
  }
else
  {
    print "Delete aborted.\n";
  }
