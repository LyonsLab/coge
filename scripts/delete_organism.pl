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

my $connstr = 'dbi:mysql:dbname=coge;host=genomevolution.org;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
my $org = $coge->resultset('Organism')->resolve($orgid);

print "Do you want to delete:  \n\t".$org->name.": ".$org->description."?\n";
my $var = <STDIN>;
if ($var =~ /y/i)
  {
    print "Initiating delete protocol.\n";
    foreach my $dsg ($org->dataset_groups)
      {
	foreach my $ds ($dsg->datasets)
	  {
	    print "Deleting dataset: ".$ds->name,"\n";
	    $ds->delete;
	  }
	my $path = $dsg->file_path;
	$path =~ s/[^\/]*$//;
	my $cmd = "rm -rf $path";
	print "Deleting genomic sequence: $cmd\n";
	`$cmd`;
	print "Deleteing dataset group\n";
	$dsg->delete;
      }
    print "Deleting organism\n";
    $org->delete;
    print "completed!\n";
  }
else
  {
    print "Delete aborted.\n";
  }
