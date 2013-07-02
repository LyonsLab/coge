#!/usr/bin/perl -w

use strict;
use CoGeX;
use DBIxProfiler;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
$coge->storage->debugobj(new DBIxProfiler());
$coge->storage->debug(1);

my $ds = $coge->resultset('Dataset')->find(6);
my $seq = $ds->get_genomic_sequence(chr=>1, start=>5000, stop=> 50000);

print length($seq),"\n";
