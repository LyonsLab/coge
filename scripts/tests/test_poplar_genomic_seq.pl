#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
$coge->storage->debug(0);

foreach my $feat ($coge->resultset('Feature')->search({'feature_names.name'=>{like=>'At1g01010'}},{join=>'feature_names'}))
  {
    print ">", join (", ", map {$_->name} $feat->feature_names),", ", $feat->type->name,": ",  $feat->start,"-", $feat->stop,'(',$feat->strand,")\n";
    print $feat->genomic_sequence,"\n";;
  }
