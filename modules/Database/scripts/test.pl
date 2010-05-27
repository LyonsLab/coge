#!/usr/bin/perl -w

use strict;
use CoGeX;
use LWP::Simple;
use Data::Dumper;


my $coge = CoGeX->dbconnect;


my $dsg = $coge->resultset('DatasetGroup')->find(7948);

my ($test) = $dsg->dataset_connectors({dataset_id=>5037});
if ($test)
  {
    print Dumper $test;
  }
		 
