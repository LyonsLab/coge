#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use CoGe::Algos::KsCalc;

my $coge = CoGeX->dbconnect();

my $fid1=shift;
my $fid2=shift;

my ($feat1) = $coge->resultset('Feature')->find($fid1);
print $feat1->genomic_sequence,"\n";
my ($feat2) = $coge->resultset('Feature')->find($fid2);
print $feat2->genomic_sequence,"\n";

my $ks = new CoGe::Algos::KsCalc();
$ks->feat1($feat1);
$ks->feat2($feat2);
my $res = $ks->KsCalc();

print Dumper $res;
