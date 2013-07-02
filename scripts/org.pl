#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
$| = 1;


my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

$s->storage->debug(1);




print $s->resultset('Dataset')->chromosomes();



