#!/usr/bin/perl -w

use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
my $GO =0;
my $dsid = shift;
unless ($dsid)
  {
    print qq{
Welcome to $0;

Usage $0 <coge database id for dataset to clone>

This program creates a new dataset entry in CoGe's database, cloning all information from the dataset specified.
};
    exit;
  }

my ($ds) = $coge->resultset('Dataset')->find($dsid);

my $new_ds = $coge->resultset('Dataset')->create({data_source_id=>$ds->data_source->id,
						 name=>$ds->name,
						 description=>$ds->description,
						 link=>$ds->link,
						 version=>$ds->version,
						 }) if $GO;
foreach my $feat ($ds->features)
  {
    print "Cloning feature: ", join ("\t", $feat->id, $feat->start, $feat->stop, $feat->strand, $feat->chromosome, $feat->feature_type_id),"\n";
    my $new_feat = $new_ds->add_to_features({feature_type_id=>$feat->feature_type->id,
					    start=>$feat->start,
					    stop=>$feat->stop,
					    strand=>$feat->strand,
					    chromosome=>$feat->chromosome,
					    }) if $GO;
    foreach my $name ($feat->feature_names)
      {
	$new_feat->add_to_feature_names({name=>$name->name,
				 description=>$name->description}) if $GO;
      }
    foreach my $loc ($feat->locations)
      {
	$new_feat->add_to_locations({start=>$loc->start,
				     stop=>$loc->stop,
				     strand=>$loc->strand,
				     chromosome=>$loc->chromosome,
				    }) if $GO;
      }
    foreach my $anno ($feat->annotations)
      {
	$new_feat->add_to_annotations({annotation=>$anno->annotation,
				       annotation_type_id=>$anno->annotation_type->id}) if $GO;
      }
    foreach my $seq ($feat->sequences)
      {
	$new_feat->add_to_sequences({sequence_type_id=>$seq->sequence_type->id,
				    sequence_data=>$seq->sequence_data}) if $GO;
      }
  }

print $ds->id, " has been cloned.\n";
print "New id for dataset object is: ",$new_ds->id,"\n" if $new_ds;
