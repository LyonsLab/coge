#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use Data::Dumper;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my ($dsid);

GetOptions (
	    "dsid=i"=>\$dsid,
	   );
unless ($dsid)
  {
    print qq{
Welcome to $0.
This program takes a dataset id from coge's genome database, uses it to find dataset_groups (with genomic sequence), and then adds those as features of type "chromosome" to the dataset.  Such "chromosome" features are necessary for coge's web-tools to work properly.

Options:  -dsid         database id for dataset

Usage:  $0 -dsid 123
};

    exit();
  }

my $ds = $coge->resultset('Dataset')->find($dsid);

unless ($ds)
  {
    print "Unable to find dataset entry for $dsid.\n";
    exit;
  }

my %chr;
foreach my $item ($ds->dataset_groups)
  {
    foreach my $gs ($item->genomic_sequences)
      {
	$chr{$gs->chromosome} = $gs->sequence_length;
      }
  }

my ($ft) = $coge->resultset('FeatureType')->search({name=>"chromosome"});
foreach my $chr (keys %chr)
  {
    my $feat = $coge->resultset('Feature')->find_or_create({
							    feature_type_id=>$ft->id,
							    dataset_id=>$ds->id,
							    start=>1,
							    stop =>$chr{$chr},
							    strand=>1,
							    chromosome=>$chr,
							    });
    my $loc = $coge->resultset('Location')->find_or_create({
							    start=>1,
							    stop =>$chr{$chr},
							    strand=>1,
							    chromosome=>$chr,
							    feature_id=>$feat->id,
							   });
    my $name = $coge->resultset('FeatureName')->find_or_create({
								name=>"chromosome $chr",
								feature_id=>$feat->id,
								});
  }
