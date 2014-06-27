#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dsid = shift;
my $feat_type_id = shift;

unless ($dsid)
  {
    print "Usage: $0 <dsid> <feature_type_id optional>\n";
    exit;
  }
my $ds = $coge->resultset('Dataset')->find($dsid);
my $ft = $coge->resultset('FeatureType')->find($feat_type_id) if $feat_type_id;
unless ($ds)
  {
    print "No dataset object for $dsid\n";
    exit;
  }

print STDERR "deleting feature for ",$ds->name,"\n";
if ($ft)
  {
    print STDERR "\tdeleting features of type ",$ft->name," only.\n";
    foreach my $feat ($ds->features({feature_type_id=>$ft->id}))
      {
	$feat->delete;
      }
  }
else
  {
    foreach my $feat ($ds->features)
      {
	    $feat->delete;
      }
  }
