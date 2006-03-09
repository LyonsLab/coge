#!/usr/bin/perl -w

use CoGe::Genome;
use strict;

my $db = new CoGe::Genome;

my $c = 0;

foreach my $fid( $db->get_feature_obj->get_all_feature_ids)
  {
    my ($feat) = $db->get_feature_obj->search(feature_id=>$fid);
    unless ($feat->locs)
      {
	$feat->delete;
	$c++;
      }
  }

print "$c features deleted\n";
