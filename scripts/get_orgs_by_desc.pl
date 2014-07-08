#!/usr/bin/perl -w

use strict;
use CoGeX;

my $desc = shift;

my $coge = CoGeX->dbconnect();

foreach my $org ($coge->resultset("Organism")->search({description=>{like=>"%".$desc."%"}}))
 {
   print join ("\t", $org->name, $org->description),"\n";
}
