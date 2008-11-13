#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;




use vars qw($DEBUG $GO $file $dsid $feat_type_name);

GetOptions (
	    "debug" => \$DEBUG,
	    "go"    => \$GO,
	    "file|f=s"=>\$file,
	    "dsid=s"=>\$dsid,
	    "feat_type_name|ftn=s" => \$feat_type_name,
	   );

unless (-r $file && $dsid && $feat_type_name)
  {
    print qq{
Usage $0 -file <tab delimited location file> -dsid <dataset database id> -feat_type_name <name of type of feature> -go

Options:
   -file | -f      name of tab delimited file with feature 

 File format:
Chromsome <tab> start <tab> stop <tab> strand <tab> name <tab> description/annotation
Names are good!  If no description/annotation is provided, none will be added.


  -feat_type_name | ftn   name of feature type e.g. gene, mRNA, CDS, CNS, gene_space

  -dsid                   database id of the dataset to which the features are added

  -go                     must be set in order to add stuff to the database

  -debug                 print messages as script is run
};
    exit();
  }


print "***************'go' flag is not set.  Nothing will be added to the database****************\n" unless $GO; 
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

print "Creating feature type $feat_type_name . . \n" if $DEBUG;
my $feat_type = $coge->resultset('FeatureType')->find_or_create({name=>$feat_type_name}) if $GO;
print "Feature type object ", $feat_type->name, " created\n" if ($GO && $feat_type && $DEBUG);

print "Creating annotation type 'annotation' . . \n" if $DEBUG;
my $anno_type = $coge->resultset('AnnotationType')->find_or_create({name=>"annotation"}) if $GO;
print "Annotation type object ", $anno_type->name, " created\n" if ($GO && $anno_type && $DEBUG);

print "Retrieving dataset for $dsid. . .\n" if $DEBUG;
my ($ds) = $coge->resultset('Dataset')->find($dsid);
    if ($ds)
  {
    print "Found dataset ", $ds->name,": ",$ds->description,"\n";
  }
else
  {
    print "Did not find a dataset for $dsid.  Please verify.\n";
    exit;
  }

#build valid chromosome list
my %chrs = map {$_=>1} $ds->chromosomes;

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
    unless ($chrs{$chr})
      {
	print "\t***$chr is not a valid chromosome for this dataset. Valid chromosomes: ";
	print join (", ", sort keys %chrs);
	print ".  Skipping this feature***\n";
	next;
      }
    print "\tCreating feature at position chr: $chr $start-$stop ($strand)\n" if $DEBUG;
    my $feat = $coge->resultset('Feature')->create(
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
