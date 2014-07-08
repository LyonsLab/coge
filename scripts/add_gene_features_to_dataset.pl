#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $GO = 0;
my $DEBUG = 1;
my $dsid;
my $add_gene =0;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

GetOptions ( "dsid=i" => \$dsid,
             "go=s"    => \$GO,
             "debug=s" => \$DEBUG,
           );

my $ds = $coge->resultset('Dataset')->find($dsid);

unless ($ds)
  {
    warn "unable to find a valid dataset entry for $dsid\n";
    exit;
  }

warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;

my ($mrna) = $coge->resultset('FeatureType')->search({name=>"mRNA"});
my ($gene) = $coge->resultset('FeatureType')->search({name=>"gene"});
foreach my $feat ($ds->features({feature_type_id=>$mrna->id}))
  {
    my $newfeat = $ds->add_to_features({
					feature_type_id => $gene->id,
					start=>$feat->start,
					stop=>$feat->stop,
					chromosome=>$feat->chromosome,
					strand=>$feat->strand,
				       }) if $GO;
    $newfeat->add_to_locations({
				start=>$feat->start,
				stop=>$feat->stop,
				chromosome=>$feat->chromosome,
				strand=>$feat->strand,
			       }) if $GO;
    foreach my $name ($feat->names)
      {
	$newfeat->add_to_feature_names({name=>$name}) if $GO;
      }
    foreach my $anno ($feat->annotations)
      {
	$newfeat->add_to_annotations({
				      annotation=>$anno->annotation,
				      annotation_type_id=>$anno->annotation_type->id}) if $GO;
	  print $anno->annotation,"::",$anno->annotation_type->id,"\n";
      }
  }
