#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my ($org, $type);
GetOptions(
	   "o|org=s"=>\$org,
	   "t|type=s"=>\$type,
	  );

unless ($org && $type)
  {
    print "Usage: $0 -o <org name> -t <feature type name>\n";
    exit;
  }

my $coge = CoGeX->dbconnect();
($org) = $coge->resultset('Organism')->resolve($org);
($type) = $coge->resultset('FeatureType')->find({name=>'CDS'});

foreach my $ds ($coge->get_current_datasets_for_org(org=>$org->id))
  {
    foreach my $feat ($ds->features({feature_type_id=>$type->id}))
      {
	print $feat->fasta,"\n";
      }

  }
