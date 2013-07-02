#!/usr/bin/perl -w
use strict;
use CoGeX;
use Data::Dumper;

my $id = shift;

unless ($id)
  {
    warn "Usage: $0 <org_id or name>\n";
    exit;
  }

my $coge = CoGeX->dbconnect();

my ($org) = $coge->resultset('Organism')->resolve($id);
foreach my $ds ($org->current_datasets())
  {
    foreach my $chr (sort $ds->chromosomes)
      {
	print join (", ", ">$chr", $org->name, $ds->name, "(v ".$ds->version.")"),"\n";
	print $ds->get_genomic_sequence(start=>1, stop=>$ds->last_chromosome_position($chr), chr=>$chr),"\n";
      }
  }
