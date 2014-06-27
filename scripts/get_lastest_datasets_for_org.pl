#!/usr/bin/perl -w

use CoGeX;
use Data::Dumper;
use strict;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $orgid = shift || 1;

foreach my $ds ($coge->get_current_datasets_for_org(org=>$orgid))
  {
    print $ds->name,"\t",$ds->version,"\n";
  }

foreach my $ds ($coge->resultset("Organism")->find($orgid)->current_datasets)
  {
    print $ds->name,"\t",$ds->version,"\n";
  }
