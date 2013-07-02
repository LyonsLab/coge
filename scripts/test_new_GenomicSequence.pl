#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use DBIxProfiler;

my $coge=CoGeX->dbconnect();
$coge->storage->debugobj(new DBIxProfiler());
$coge->storage->debug(1);

foreach my $ds ($coge->resultset('Dataset')->all())
  {
    foreach my $chr ($ds->chromosomes())
      {
	print $ds->name,"\t",$ds->get_genomic_sequence(chr=>$chr),"\n";
      }
  }
