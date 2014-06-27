#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my $dsid1 = shift;
my $dsid2 = shift;

unless ($dsid1 && $dsid2)
  {
    print qq{
Usage:
$0 <dataset id 1> <dataset id 2>

This program will copy the features from dataset 1 to dataset 2.
};
    exit;
  }

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my ($ds1) = $coge->resultset('Dataset')->resolve($dsid1);
my ($ds2) = $coge->resultset('Dataset')->resolve($dsid2);

unless ($ds1 && $ds2)
  {
    print "Problem with retrieving one of the datasets.\n  Please check their names or ids.  Exiting. . .";
    exit;
  }

foreach my $feat1 ($ds1->features)
  {
    my $feat2 = $ds2->add_to_features(
#    print Dumper
				      {feature_type_id=>$feat1->type->id,
				       start=>$feat1->start,
				       stop=>$feat1->stop,
				       strand=>$feat1->strand,
				       chromosome=>$feat1->chromosome,
#				      };
				      });
    foreach my $name($feat1->feature_names)
      {
	$feat2->add_to_feature_names(
#	print Dumper
				     {
				      name=>$name->name,
				      description=>$name->description,
				      primary_name=>$name->primary_name,
#				     };
				     });
      }
    foreach my $loc($feat1->locations)
      {
	$feat2->add_to_locations(
#	print Dumper
				 {
				  start=>$loc->start,
				  stop=>$loc->stop,
				  strand=>$loc->strand,
				  chromosome=>$loc->chromosome,
#				 };
				});
      }
    foreach my $anno($feat1->annotations)
      {
	$feat2->add_to_annotations(
#	print Dumper
				   {
				    annotation=>$anno->annotation,
				    annotation_type_id=>$anno->annotation_type_id,
#				   };
				   });
      }
    foreach my $seq($feat1->sequences)
      {
	$feat2->add_to_sequences (
#	print Dumper
				  {
				   sequence_type_id=>$seq->sequence_type_id,
				   sequence_data=>$seq->sequence_data,
#				  };
				  });

      }
  }
