#!/usr/bin/perl

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=localhost;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $file = "/home/elyons/projects/CoGeX/tmp/new_genomic_sequence_table_data.txt";

open (IN, $file);
while (<IN>)
  {
    chomp;
    next if /^#/;
    next unless $_;
    my ($ds_id, $start, $stop, $chr, $seq_type_id) = split/\t/;
    unless ($ds_id && $start && $stop && $chr && $seq_type_id) {warn "Missing data: $_\n"; next;}
#    unless ($stop > 1) {warn "Stop problem: $_\n";}
    $coge->resultset('GenomicSequence')->create({start=>$start,
						 stop=>$stop,
						 chromosome=>$chr,
						 dataset_id=>$ds_id,
						 genomic_sequence_type_id=>$seq_type_id});
  }
close IN;
