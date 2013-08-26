#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $anno_type = shift;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my %data;
foreach my $anno ($coge->resultset('Annotation')->search({annotation_type_id=>$anno_type}))
  {
    if (    $data{$anno->annotation_type_id}{$anno->feature_id}{$anno->annotation})
      {
	print "Duplicate annotation found\n";
	$anno->delete;
	next;
      }
    $data{$anno->annotation_type_id}{$anno->feature_id}{$anno->annotation}=1;
  }
