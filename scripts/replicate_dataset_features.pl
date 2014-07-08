#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;

my ($dsid1, $dsid2, $ftid, @skip_ftids, $db, $user, $pass);

GetOptions(
	   "dsid1=i"=>\$dsid1,
	   "dsid2=i"=>\$dsid2,
	   "ftid=i"=>\$ftid,
	   "skip_ftid=i"=>\@skip_ftids,
             "database|db=s"=>\$db,
            "user|u=s"=>\$user,
            "password|pw=s"=>\$pass,

	   );

unless ($dsid1 && $dsid2)
  {
    print qq{
Usage:
$0 -dsid1 <dataset id 1> -dsid2 <dataset id 2>

This program will copy the features from dataset 1 to dataset 2.

Options:

 -dsid1              database dataset id for the dataset from which features are copied

 -dsid2              database dataset id for the dataset to which features are copied

 -ftid               OPTIONAL:  feature type id for the features to be copied. This is useful if, for example,
                                only "chromosome" (ftid 4) features are to be copied.  If left undefined, all
                                features will be copied.

-skip_ftid          OPTION:  skip features of a particular type.  E.g. chromosomes

 -db                coge database name

 -user              database user

 -pw                database password

};
    exit;
  }

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3306";
my $coge = CoGeX->connect($connstr, $user, $pass );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my ($ds1) = $coge->resultset('Dataset')->find($dsid1);
my ($ds2) = $coge->resultset('Dataset')->find($dsid2);

unless ($ds1 && $ds2)
  {
    print "Problem with retrieving one of the datasets.\n  Please check their names or ids.  Exiting. . .";
    exit;
  }

my $search = {feature_type_id=>$ftid} if $ftid;
my %skip_ids = map{$_=>1} @skip_ftids;

my $count =0;
my %types;
print "Starting replication from ",$ds1->dataset_groups->[0]->organism->name,": ",$ds1->name, " to ", $ds1->dataset_groups->[0]->organism->name, ": ",$ds2->name, "\n";
foreach my $feat1 ($ds1->features($search))
  {
    next if $skip_ids{$feat1->type->id};
    print ".";
    print "\n" unless $count % 100;
    $types{$feat1->type->id}++;
#    sleep (.3);
    my $chr = $feat1->chromosome;
    my $feat2 = $ds2->add_to_features(
#    print Dumper
				      {feature_type_id=>$feat1->type->id,
				       start=>$feat1->start,
				       stop=>$feat1->stop,
				       strand=>$feat1->strand,
				       chromosome=>$chr,
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
#				  chromosome=>$loc->chromosome,
				  chromosome=>$chr,
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
    $count++;
  }
print $count, "Features Processed\n";
print Dumper \%types;
