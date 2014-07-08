#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
$| = 1;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

$s->storage->debug(1);

print $s->resultset('Dataset')->chromosomes();
