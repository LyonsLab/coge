#!/usr/bin/perl -w
use strict;
use CoGeX;


my $source_id = shift || 21;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my %data;

my $source = $coge->resultset('DataSource')->find($source_id);
foreach my $ds ($source->datasets)
  {
    if ($data{uc($ds->name)})
      {
	print $ds->name ," exists!\n";
	$ds->delete;
      }
    else
      {
	$data{uc($ds->name)}=1;
      }
  }

