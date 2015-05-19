#!/usr/bin/perl -w

use CoGeX;
use strict;

my $connstr = 'dbi:mysql:dbname=coge;host=genomevolution.org;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
my $GO =1;

my $dsgid = shift;

my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);

my $new_dsg = $coge->resultset('DatasetGroup')->create({
							name=>$dsg->name." CLONE",
							description=>$dsg->description,
							version=>$dsg->version,
							organism_id=>$dsg->organism_id,
							genomic_sequence_type_id=>$dsg->genomic_sequence_type_id,
							file_path=>$dsg->file_path}) if $GO;
print "Cloned dataset group id: ".$new_dsg->id,"\n" if $new_dsg;
foreach my $gs ($dsg->genomic_sequences)
  {
#    print "cloning chromosome: ".$gs->chromosome,"\n";
    $new_dsg->add_to_genomic_sequences({
					sequence_length=>$gs->sequence_length,
					chromosome=>$gs->chromosome,
					}) if $GO;
  }

foreach my $ds ($dsg->datasets)
  {
    my $new_ds = $coge->resultset('Dataset')->create({
						      name=>$ds->name. " CLONE",
						      description=>$ds->description,
						      version=>"clone",
						      data_source_id => $ds->data_source_id,
						      }) if $GO;
    print "Cloned dataset id: ".$new_ds->id,"\n" if $new_ds;
    $coge->resultset('DatasetConnector')->create({
						  dataset_id=>$new_ds->id,
						  dataset_group_id=>$new_dsg->id,
						  }) if $GO;
    foreach my $chr ($ds->features({feature_type_id=>301}) )#get the chromosomes
      {
	my $feat = $new_ds->add_to_features({start=>$chr->start,
					     stop=>$chr->stop,
					     strand=>$chr->strand,
					     chromosome=>$chr->chromosome,
					     feature_type_id=>$chr->feature_type_id,
					    }) if $GO;
  	foreach my $fn ($chr->feature_names)
	  {
	    $feat->add_to_feature_names({name=>$fn->name,
					 description=>$fn->description,
					}) if $GO;
	  }
  	foreach my $loc ($chr->locations)
	  {
	    $feat->add_to_locations({start=>$loc->start,
				     stop=>$loc->stop,
				     strand=>$loc->strand,
				     chromosome=>$loc->chromosome,
				    }) if $GO;
	  }
	foreach my $anno ($chr->annotations)
	  {
	    $feat->add_to_annotations({annotation=>$anno->annotation,
				       annotation_type_id=>$anno->annotation_type_id,
				      }) if $GO;
	  }
      }
  }
