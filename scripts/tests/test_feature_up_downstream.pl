#!/usr/bin/perl -w

use CoGeX;
use Data::Dumper;

my $fid = 25 || shift;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $feat = $coge->resultset('Feature')->find($fid);

print join ("\t", $feat->start, $feat->stop, $feat->strand, $feat->chromosome, $feat->dataset->id, $feat->dataset->name),"\n";

print join ("\n", map {"\t".join ("\t", $_->start,$_->stop, $_->strand,$_->chromosome)} $feat->locations),"\n";
print $feat->genomic_sequence,"\n";
print $feat->genomic_sequence(up=>10,down=>10),"\n";
