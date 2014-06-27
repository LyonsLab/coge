#!/usr/bin/perl -w

use CoGe::Genome;
use Data::Dumper;
use strict;

my $db = new CoGe::Genome;

unless (@ARGV)
  {
    print qq{
USAGE: $0 <data_information id from the CoGe Genome database>
};
    exit;
  }
foreach my $id (@ARGV)
  {
    my ($di) = $db->get_data_info_obj->search({data_information_id=>$id});
    unless ((ref $di) =~ /CoGe/i)
      {
	print "$id didn't return a CoGe object . . . skipping\n";
	next;
      }
    print "Deleting data_information: ",$di->name,": ", $di->description, " version ", $di->version,". . .";
    $di->delete;
    print "Done\n";
  }
