#!/usr/bin/perl

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

print "#";
print
  join( "\t", qw{dataset_id start stop chromosome genomic_sequence_type_id} ),
  "\n";
foreach my $item (
    $coge->resultset('Dataset')->search( {}, { order_by => 'dataset_id' } ) )
{
    next unless $item->sequence_type;
    my $dsid      = $item->id;
    my $start     = 1;
    my $type_id   = $item->sequence_type->id;
    my $chr_count = $item->chromosomes;
    foreach my $chr ( $item->chromosomes ) {
        my $stop = $item->last_chromosome_position($chr);
        print join( "\t", $dsid, $start, $stop, $chr, $type_id ), "\n";
    }
}
