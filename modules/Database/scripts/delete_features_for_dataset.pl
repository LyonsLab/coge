#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);


my $dsid = shift;
unless ($dsid)
  {
    print "Usage: $0 <dsid>\n";
    exit;
  }
my $ds = $coge->resultset('Dataset')->find($dsid);
unless ($ds)
  {
    print "No dataset object for $dsid\n";
    exit;
  }

print STDERR "deleting feature for ",$ds->name,"\n";
foreach my $feat ($ds->features)
  {
    $feat->delete;
  }

