#!/usr/bin/perl -w
use strict;
use CoGeX;
use Data::Dumper;

my $id = shift;

unless ($id)
  {
    warn "Usage: $0 <dataset_id>\n";
    exit;
  }

my $coge = CoGeX->dbconnect();
my $res = $coge->resultset('Dataset')->search({dataset_id=>$id});
my $ds = $res->next();

foreach my $chr (sort $ds->chromosomes) {
    print ">$chr\n";
    print $ds->get_genomic_sequence(start=>1, stop=>$ds->last_chromosome_position($chr), chr=>$chr) . "\n";
}
