#!/usr/bin/perl -w
use strict;
use CoGeX;
use Data::Dumper;
use DBIxProfiler;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);


my $org = shift;
($org) = $coge->resultset('Organism')->resolve($org);
print join ("\n", map {$_->name} $org->genomic_sequence_types),"\n";
my ($type) = $coge->resultset('GenomicSequenceType')->resolve(2);
foreach my $ds ($org->current_datasets(type=>$type))
  {
    print $ds->name,": ",$ds->sequence_type->name,"\n";
#    my ($gs) = $ds->genomic_sequences;
#    print $gs->genomic_sequence_type->name,"\n";
  }
