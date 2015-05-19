#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

use vars qw($dsid $GO $DEBUG);

GetOptions ( "dsid=i" => \$dsid,
             "go=s"    => \$GO,
             "debug=s" => \$DEBUG,
           );

unless ($dsid)
  {
    print qq{
Welcome to $0.

Usage:  $0 -dsid <dataset id for annotations> -go 1 < file_with_annotations

Use to go flag to run.  Keeps you from running accidentally.

};
    exit;
  }

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $ds = $coge->resultset('Dataset')->find($dsid);

my $feat_type = $coge->resultset('FeatureType')->find_or_create(name=>"manually predicted pre-tRNA", description=>"Blastn assisted");

my $anticodon_type = $coge->resultset('AnnotationType')->find_or_create(name=>"anticodon");
my $codon_type =  $coge->resultset('AnnotationType')->find_or_create(name=>"codon");

while (<>)
  {
    chomp;
    my @line = split /\t/;
    my $start = $line[1];
    my $stop = $start+$line[2]-1;
    my $strand = 1;
    my $feat = $ds->add_to_features({start=>$start,
				    stop=>$stop,
				    chromosome=>1,
				    strand=>$strand,
				    feature_type_id=>$feat_type->id}) if $GO;
    my $loc = $feat->add_to_locations({start=>$start,
				      stop=>$stop,
				      chromosome=>1,
				      strand=>$strand,
				      }) if $GO;
    my $anno = $feat->add_to_annotations({annotation=>"Tyrosine, Tyr GTA", annotation_type_id => $anticodon_type->id}) if $GO;
    $anno = $feat->add_to_annotations({annotation=>"Tyrosine, Tyr TAC", annotation_type_id => $codon_type->id}) if $GO;
    $feat->add_to_feature_names({name=>$line[0]."_oa", description=>"reannotation of incomplete gene model"}) if $GO;
  }
