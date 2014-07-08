#!/usr/bin/perl -w

use strict;
use CoGeX;
use LWP::Simple;
use Data::Dumper;

my $coge = CoGeX->dbconnect;

foreach my $org ($coge->resultset('Organism')->all)
  {
    print $org->name,"\n";
  }
