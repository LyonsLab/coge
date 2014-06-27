#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my ($org, $type);
GetOptions(
	   "o|org=s"=>\$org,
	   "t|type=s"=>\$type,
	  );
$type = "CDS";
unless ($org && $type)
  {
    print "Usage: $0 -o <org name> -t <feature type name>\n";
    exit;
  }

my $coge = CoGeX->dbconnect();
($org) = $coge->resultset('Organism')->resolve($org);
($type) = $coge->resultset('FeatureType')->find({name=>'CDS'});

my @data;
foreach my $ds ($coge->get_current_datasets_for_org(org=>$org->id))
  {
    foreach my $feat ($ds->features({feature_type_id=>$type->id}))
      {
	my ($name) = $feat->names;
#	print join ("\t", $name, map {sprintf("%.2f",100*$_)} $feat->wobble_content),"\n";
	push @data, [$name, map {sprintf("%.2f",100*$_)} $feat->wobble_content];
      }

  }
print join ("\t", qw(NAME GC AT)),"\n";
my $limit =100;
my $count = 0;
foreach my $item (sort {$a->[1]<=>$b->[1]} @data)
  {
    print join ("\t", @$item),"\n";
    $count++;
    last if $count == 100;
  }
print "\n\n";
$count = 0;
foreach my $item (sort {$b->[1]<=>$a->[1]} @data)
  {
    print join ("\t", @$item),"\n";
    $count++;
    last if $count == 100;
  }
