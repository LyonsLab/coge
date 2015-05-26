#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

use vars qw($DEBUG $GO $file $dsid $feat_type_name  $ds_name $ds_desc $ds_link $ds_version $source_name $source_desc $source_link $source_id);

GetOptions (
	    "source_name=s" => \$source_name, # datasource
	    "source_desc=s" => \$source_desc,
	    "source_link=s" => \$source_link,
	    "source_id=s"   => \$source_id,
	    "ds_name=s" => \$ds_name,# datasetid
	    "ds_desc=s" => \$ds_desc,
	    "ds_link=s" => \$ds_link,
	    "ds_version=s" => \$ds_version,
	    "dsid=i" => \$dsid,
	    "debug" => \$DEBUG,
	    "go"    => \$GO,
	    "file|f=s"=>\$file,
	    "feat_type_name|ftn=s" => \$feat_type_name,
	   );

unless ($file && -r $file && $feat_type_name)
  {
    usage();
    exit();
  }

print "***************'go' flag is not set.  Nothing will be added to the database****************\n" unless $GO;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

print "Creating feature type $feat_type_name . . \n" if $DEBUG;
my $feat_type = $coge->resultset('FeatureType')->find_or_create({name=>$feat_type_name}) if $GO;
print "Feature type object ", $feat_type->name, " created\n" if ($GO && $feat_type && $DEBUG);

print "Creating annotation type 'annotation' . . \n" if $DEBUG;
my $anno_type = $coge->resultset('AnnotationType')->find_or_create({name=>"annotation"}) if $GO;
print "Annotation type object ", $anno_type->name, " created\n" if ($GO && $anno_type && $DEBUG);

print "Creating source for dataset (if necessary)\n";
if ($source_name)
  {
    my $source = $coge->resultset("DataSource")->find_or_create({name=>$source_name,description=>$source_desc, link=>$source_link});
    $source_id = $source->id;
  }

print "Creating or retrieving dataset\n";
my $ds = generate_ds(ds_name => $ds_name,
		     ds_desc => $ds_desc,
		     ds_link => $ds_link,
		     ds_version => $ds_version,
		     ds_id =>$dsid,
		     source_id=>$source_id,
		    );
print "Retrieving dataset for $dsid. . .\n" if $DEBUG;
if ($ds)
  {
    print "Found dataset ", $ds->name,": ",$ds->description,"\n";
  }
else
  {
    print "Did not find a dataset for $dsid.  Please verify.\n";
    usage();
    exit;
  }

#build valid chromosome list
#my %chrs = map {$_=>1} $ds->chromosomes;

open (IN, $file);
while (<IN>)
  {
    next if /^#/;
    chomp;
    next unless $_;
    my ($chr, $start, $stop, $strand, $name, $desc) = split/\t/;
#    s/;/_/g;
    if ($DEBUG)
      {
	print "Working on: ";
	print join ("\t", $chr, $start, $stop, $strand, $name);
	print "\t", $desc if $desc;
	print "\n";
      }
#     unless ($chrs{$chr})
#       {
# 	print "\t***$chr is not a valid chromosome for this dataset. Valid chromosomes: ";
# 	print join (", ", sort keys %chrs);
# 	print ".  Skipping this feature***\n";
# 	next;
#       }
    print "\tCreating feature at position chr: $chr $start-$stop ($strand)\n" if $DEBUG;
    my $feat = $coge->resultset('Feature')->find_or_create(
							   {feature_type_id=>$feat_type->id,
							    dataset_id=>$ds->id,
							    start=>$start,
							    stop=>$stop,
							    strand=>$strand,
							    chromosome=>$chr,
							   }) if $GO;
    print "\t\tCreating location chr: $chr $start-$stop ($strand)\n" if $DEBUG;
    $feat->add_to_locations({
			     start=>$start,
			     stop=>$stop,
			     strand=>$strand,
			     chromosome=>$chr,
			     }) if $GO;

    print "\t\tAdding name to feature: $name\n" if $DEBUG;
    $feat->add_to_feature_names({
				 name=>$name,
				 }) if $GO;
    print "\t\tAdding annotation to feature: $desc\n" if $DEBUG && $desc;
    $feat->add_to_annotations({annotation=>"$desc",
			       annotation_type_id=>$anno_type->id,
			       }) if $desc && $GO;
  }
close IN;

print "Completed loading feature to dataset: ".$ds->name. " (database id: ".$ds->id.")\n";

sub generate_ds
  {
    my %opts = @_;
    my $ds_name = $opts{ds_name};
    my $ds_desc = $opts{ds_desc};
    my $ds_link = $opts{ds_link};
    my $ds_version = $opts{ds_version};
    my $ds_id = $opts{ds_id};
    my $source_id = $opts{source_id};
    unless ($ds_name || $ds_id)
      {
	warn "no dataset name or database id specified\n";
	return;
      }
    my $ds = $ds_id ? $coge->resultset('Dataset')->find($ds_id) :
      $coge->resultset('Dataset')->find_or_create({
						   name                => $ds_name,
						   description         => $ds_desc,
						   link                => $ds_link,
						   data_source_id      => $source_id,
						   version=>$ds_version,
						  });;
    return $ds;

  }

sub usage
  {

    print qq{
Usage $0 -file <tab delimited location file> -dsid <dataset database id> -feat_type_name <name of type of feature> -go

Options:
   -file | -f      name of tab delimited file with feature

 File format:
Chromsome <tab> start <tab> stop <tab> strand <tab> name <tab> description/annotation
Names are good!  If no description/annotation is provided, none will be added.

  -feat_type_name | ftn   name of feature type e.g. gene, mRNA, CDS, CNS, gene_space

  -go                     must be set in order to add stuff to the database

  -debug                 print messages as script is run

  -dsid                   database id of the dataset to which the features are added

########################################
#  IF YOU NEED TO CREATE A DATASET:
#
#  You will need to either create a "source" for your dataset, or specify the database id of the source
########################################

 -source_name        name of source of data (e.g. "Freeling Lab")

 -source_desc        description of source of data (e.g. "Best lab in the world")

 -source_link        html/url link to the source (e.g. "http://microscopy.berkeley.edu/~freeling/")

 -source_id          database id of source of data (e.g. "69")

####################
#  Create a dataset
####################
 -ds_name           name of dataset

 -ds_desc           description of dataset

 -ds_link           html/url link to dataset

 -ds_version        version of dataset

};

  }
