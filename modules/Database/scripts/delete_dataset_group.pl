#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my $connstr = 'dbi:mysql:dbname=coge;host=synteny.cnr.berkeley.edu;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
my ($dsgid, $delete_seqs);

GetOptions (
	    "dsgid=i"=>\$dsgid,
	    "delete_seqs"=>\$delete_seqs
	    );


unless ($dsgid)
  {
    print qq{
welcome to $0

Usage:  $0 -dsgid <database id for dataset group> -delete_seqs

This will delete the database entry for the dataset group, associated entries in the database_connector table, and the datasets if the datasets do not belong to another dataset group!  

To delete the genomic sequences associted with the dataset group, add flag -delete_seqs
  };
    exit;
  }


my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
unless ($dsg)
  {
    print "unable to find entry for $dsgid\n";
    exit;
  }

delete_dsg($dsg);
print "Finished deleting dataset group\n";

sub delete_dsg
  {
    my $dsg = shift;
    foreach my $ds ($dsg->datasets)
      {
	my %dsgs;
	foreach my $item ($ds->dataset_groups)
	  {
	    next if $item->id eq $dsg->id;
	    $dsgs{$item->id}=$item;
	  }
	next if keys %dsgs;
	print "Deleting dataset ".$ds->name," (".$ds->id.")"."\n";
	$ds->delete();
      }
    if ($delete_seqs)
      {
	my $path = $dsg->file_path;
	$path =~ s/[^\/]*$//;
	my $cmd = "rm -rf $path";
	print "Deleting genomic sequence: $cmd\n";
	`$cmd`;
      }
    else
      {
	print "DID NOT DELETE GENOMIC SEQUENCES: ". $dsg->file_path,"\n";
      }
    $dsg->delete();
  }

