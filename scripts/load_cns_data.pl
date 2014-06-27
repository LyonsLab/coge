#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
my $GO = 1; #set to 1 to load into database

my $file = "/tmp/CNSwithCompleteInfo";

my $org = $coge->resultset('Organism')->resolve(1);
my %ds;
my $version;
foreach my $ds ($org->current_datasets)
  {
    my ($chr) = $ds->get_chromosomes();
    $ds{$chr}=$ds;
    $version = $ds->version unless $version;
    if ($version && $version ne $ds->version)
      {
        warn "Problem with dataset version:  previously identified: $version, new: ",$ds->version,"\n";
      }
  }

open (IN, $file);
#<IN>; #toss header
my $feat_type = $coge->resultset('FeatureType')->find_or_create({name=>"CNS"}) if $GO;
my $anno_type = $coge->resultset('AnnotationType')->find_or_create({name=>"annotation"}) if $GO;

while (<IN>)
  {
    chomp;
    next unless $_;
    my ($name, $gene, $strand, $chr, $start, $stop) = split/\t/;
#    s/;/_/g;
    print join ("\t", $name, $gene, $strand, $chr, $start, $stop),"\n";;
    my $feat = $coge->resultset('Feature')->create(
						   {feature_type_id=>$feat_type->id,
						    dataset_id=>$ds{$chr}->id,
						    start=>$start,
						    stop=>$stop,
						    strand=>$strand,
						    chromosome=>$chr,
						   }) if $GO;
    $feat->add_to_locations({
			     start=>$start,
			     stop=>$stop,
			     strand=>$strand,
			     chromosome=>$chr,
			     }) if $GO;
    $feat->add_to_feature_names({
				 name=>$name,
				 }) if $GO;
    $feat->add_to_feature_names({
				 name=>$gene."_cns",
				 }) if $GO;
    $feat->add_to_annotations({annotation=>"CNS for gene $gene",
			       annotation_type_id=>$anno_type->id,
			       }) if $GO;
  }
close IN;
