#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my $connstr = 'dbi:mysql:dbname=coge;host=genomevolution.org;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my ($dsid, $file, $GO);
GetOptions (
            "dsid=i" => \$dsid,
             "go=s"    => \$GO,
	   "file=s"=>\$file,
           );

unless ($dsid)
  {
    print qq{
Welcome to $0!

Usage:  $0 -dsid <coge dataset id to which annotations are added> -file <tRNA-scan tab-delimited output> -go 1

 without -go 1, no data will be loaded

};
    exit;
  }

my $ds = $coge->resultset('Dataset')->find($dsid);
my ($tRNA_feat_type) = $coge->resultset('FeatureType')->search({name=>'tRNA'});
my ($gene_feat_type) = $coge->resultset('FeatureType')->search({name=>'gene'});
my ($pseudogene_feat_type) = $coge->resultset('FeatureType')->search({name=>'pseudogene'});
my $product_anno_type = $coge->resultset('AnnotationType')->find_or_create({name=>"product"});
my $anticodon_anno_type = $coge->resultset('AnnotationType')->find_or_create({name=>"Anticodon"});
my $note_anno_type = $coge->resultset('AnnotationType')->find_or_create({name=>"Note"});
open (IN, $file) || die "can;t open $file for reading";
while (<IN>)
  {
    next if /^Sequence/;
    next if /^Name/;
    next if /^-----/;
    chomp;
    my @line = split /\s+/;
    my $strand = $line[2] > $line[3] ? -1 : 1;
    my ($start, $stop) = ($line[2], $line[3]);
    ($start, $stop) = ($stop, $start) if $start > $stop;
    my @feats;
    if ($line[4] eq "Pseudo")
      {
	push @feats, $ds->add_to_features({
					   feature_type_id => $pseudogene_feat_type->id,
					   start=>$start,
					   stop=>$stop,
					   chromosome=>1,
					   strand=>$strand,
					  }) if $GO;
      }
    else
      {
	push @feats, $ds->add_to_features({
					   feature_type_id => $gene_feat_type->id,
					   start=>$start,
					   stop=>$stop,
					   chromosome=>1,
					   strand=>$strand,
					  }) if $GO;
	push @feats, $ds->add_to_features({
					   feature_type_id => $tRNA_feat_type->id,
					   start=>$start,
					   stop=>$stop,
					   chromosome=>1,
					   strand=>$strand,
					  }) if $GO;
      }
    foreach my $feat (@feats)
#    my $feat;
      {
	print "working on feature ".$feat->id,"\n";
	$feat->add_to_locations(
				{
				 start      => $start,
				 stop       => $stop,
				 strand     => $strand,
				 chromosome => 1
				}
			       ) if $GO;
	my $name = "tRNA-".$line[4];
#	print $name,"\n";
	$feat->add_to_feature_names({
				 name=>$name,
				}) if $GO;

	$feat->add_to_annotations({annotation=>$name, annotation_type_id => $product_anno_type->id}) if $GO;
	my $anticodon = $line[5];
	$anticodon =~ s/T/U/g;
#	print $anticodon,"\n";
	$feat->add_to_annotations({annotation=>$line[5], annotation_type_id => $anticodon_anno_type->id}) if $GO;
	$feat->add_to_annotations({annotation=>"Predicted by tRNAscan-SE v1.23",link=>"ftp://selab.janelia.org/pub/software/tRNAscan-SE/", annotation_type_id => $note_anno_type->id}) if $GO;
      }
  }
close IN;
