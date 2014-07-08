#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;
use File::Path;

my ($dsgid, $new_seq_dir, $db1, $db2, $u1, $u2, $p1, $p2);

#./replicate_genome_between_coge_installations.pl -dsgid 1149 -sd /opt/apache/Oryza_CoGe/data/genomic_sequence/
GetOptions(
	   "dsgid=i"=>\$dsgid,
	   "sequence_dir|sd=s"=>\$new_seq_dir,
	   "u1=s"=>\$u1,
	   "u2=s"=>\$u2,
	   "p1=s"=>\$p1,
	   "p2=s"=>\$p2,
	   "db1=s"=>\$db1,
	   "db2=s"=>\$db2,
	   );

$u2 = $u1 unless $u2;
$p2 = $p1 unless $p2;
$db2 = $db1 unless $db2;

unless ($dsgid && $u1 && $u2 && $db1 && $db2 && $p1 && $p2 && $new_seq_dir)
  {
    help();
  }

my $connstr1 = "dbi:mysql:dbname=$db1;host=localhost;port=3306";
my $connstr2 = "dbi:mysql:dbname=$db2;host=localhost;port=3306";
my $coge1 = CoGeX->connect($connstr1, $u1, $p1 );
my $coge2 = CoGeX->connect($connstr2, $u2, $p2 );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dsg = $coge1->resultset('DatasetGroup')->find($dsgid);

unless ($dsg)
  {
    print "Problem with retrieving the coge database object for the genometo be copied.\nExiting.\n";
    exit;
  }

print "Finding or creating organism object\n";
my $org = $coge2->resultset('Organism')->find_or_create({name=>$dsg->organism->name,
							 description=>$dsg->organism->description});
print "Finding or creating genomic sequence type object\n";
my $gst = $coge2->resultset('GenomicSequenceType')->find_or_create({name=>$dsg->type->name,
								    description=>$dsg->type->description});
print "Replicating dataset_group_object\n";
my $new_dsg = $org->add_to_dataset_groups({name=>$dsg->name,
					   description=>$dsg->desc,
					   version=>$dsg->version,
					   organism_id => $org->id,
					   genomic_sequence_type_id=>$gst->id,
					  });
replicate_sequences(dsg1=>$dsg, dsg2=>$new_dsg);

foreach my $gs ($dsg->genomic_sequences)
  {
    $new_dsg->add_to_genomic_sequences({
					sequence_length=>$gs->sequence_length,
					chromosome=>$gs->chromosome,
				       });
  }
my @dss;
print "Replicating datasets\n";
foreach my $ds ($dsg->datasets)
  {
    my $source = $coge2->resultset('DataSource')->find_or_create({
								  name=>$ds->source->name,
								  description=>$ds->source->desc,
								  link=>$ds->source->link,
								 });
    my $new_ds = $coge2->resultset('Dataset')->create({
						      name=>$ds->name,
						      description => $ds->description,
						      version=>$ds->version,
						      link=>$ds->link,
						      data_source_id=>$source->id,
						     });
    #create linker
    $coge2->resultset('DatasetConnector')->create({
						  dataset_id=>$new_ds->id,
						  dataset_group_id=>$new_dsg->id,
						 });
    push @dss, [$ds, $new_ds];
  }
print "Replicating features\n";
foreach my $item (@dss)
  {
    my ($ds1, $ds2) = @$item;

    foreach my $feat1 ($ds1->features)
      {
	my $ft = $coge2->resultset('FeatureType')->find_or_create({name=>$feat1->type->name,
								   description=>$feat1->type->desc});
	my $feat2 = $ds2->add_to_features(
					  #    print Dumper
					  {feature_type_id=>$ft->id,
					   start=>$feat1->start,
					   stop=>$feat1->stop,
					   strand=>$feat1->strand,
					   chromosome=>$feat1->chromosome,
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
				     });
	  }
	foreach my $anno($feat1->annotations)
	  {
	    my $at_search = {name=>$anno->type->name,
			     description=>$anno->type->desc};
	    if ($anno->type->group)
	      {
		my $atg = $coge2->resultset('AnnotationTypeGroup')->find_or_create({name=>$anno->type->group->name,
										    description=>$anno->type->group->desc});
		$at_search->{annotation_type_group_id}=$atg->id;
	      }
	    my $at = $coge2->resultset("AnnotationType")->find_or_create($at_search);
	    $feat2->add_to_annotations(
				       #	print Dumper
				       {
					annotation=>$anno->annotation,
					annotation_type_id=>$at->id,
					link=>$anno->link,
					#				   };
				       });
	  }
	foreach my $seq($feat1->sequences)
	  {
	    my $st = $coge2->resultset('SequenceType')->find_or_create({name=>$seq->type->name,
									description=>$seq->type->desc});
	    $feat2->add_to_sequences (
				      #	print Dumper
				      {
				       sequence_type_id=>$st->id,
				       sequence_data=>$seq->sequence_data,
				       #				  };
				      });
	  }
      }
  }

print "Finished\n";
print "New DatasetGroup id: ".$new_dsg->id,"\n";
print "Link:\n";
print "\t"."http://genomevolution.org/CoGe/GenomeView.pl?dsgid=".$new_dsg->id."\n";

sub replicate_sequences
  {
    my %opts = @_;
    my $dsg1 = $opts{dsg1};
    my $dsg2 = $opts{dsg2};

    print "Replicating sequence data\n";
    my $seq_path = $dsg1->file_path;
    $seq_path =~ s/[^\/]*$//;
    my $new_path = $dsg2->get_path;
    $new_path = $new_seq_dir."/".$new_path;
    mkpath $new_path;
    print "\tOld path: ".$seq_path,"\n";
    print "\tNew path: ".$new_path,"\n";
    `cp -r $seq_path/* $new_path`;

    my $dsgid1 = $dsg1->id;
    my $dsgid2 = $dsg2->id;
    opendir (DIR, $new_path);
    while (my $item = readdir(DIR))
      {
	next unless $item =~ /$dsgid1/;
	print $item,"\n";
	my $new_name = $item;
	$new_name =~ s/$dsgid1/$dsgid2/;
	my $cmd = "mv $new_path/$item $new_path/$new_name";
	print "\t".$cmd,"\n";
	`$cmd`;
      }
    closedir DIR;
    $new_path .= "/".$dsg2->id.".faa";
    $dsg2->file_path($new_path);
    $dsg2->update;

  }

sub help
  {
    print qq{
Usage:
$0 -dsgid <dataset_group id>

This program generates a copy of a coge dataset_group (genome).  The genomic sequence is NOT copied and that from the original genome is used.

Options:

 -dsgid              coge database dataset_group

 -sequence_dir | sd  directory to which sequence directory is replicated

 -db1               source coge database name
 -u1                source coge database user name
 -p1                source coge database user password

 -db2               sink coge database name
 -u2                sink coge database user name
 -p2                sink coge database user password

If not specified, db1, u1, and p1 will be used for the sink.

Example:

 ./replicate_genome_between_coge_installations.pl -sd /opt/apache/Oryza_CoGe/data/genomic_sequence/ -dsgid 1149

};
    exit;
  }
