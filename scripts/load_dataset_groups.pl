#!/usr/bin/perl

use strict;
use Data::Dumper;
use DBI;

my $file = '/home/elyons/projects/CoGeX/tmp/dataset_groups.txt';
my $connstr = 'dbi:mysql:dbname=DB;host=localhost;port=PORT';
my $dbi = DBI->connect($connstr, 'USER', 'PASSWORD' );
my $group_insert = $dbi->prepare(qq{
INSERT INTO dataset_group (version, organism_id, genomic_sequence_type_id) VALUES (?,?,?)
});
my $connector_insert = $dbi->prepare(qq{
INSERT INTO dataset_connector (dataset_id, dataset_group_id) VALUES (?,?)
});

process_file($file);
sub process_file
  {
    my $file = shift;
    open (IN, $file);
    while (<IN>)
      {
	next if /^#/;
	chomp;
	next unless $_;
	my ($ver, $org_id, $seq_type_id, $datasets) = split/\t/;
	my @datasets = split /::/, $datasets;
	my $val = $group_insert->execute($ver, $org_id, $seq_type_id);
	my $dataset_group_id = $dbi->last_insert_id(undef, undef, undef, undef);
	foreach my $dsid (@datasets)
	  {
	    $connector_insert->execute($dsid, $dataset_group_id);
	  }
      }
  }
