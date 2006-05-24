#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGe::Genome;
use LWP::Simple;
use POSIX;

my ($org_id, $version, $zoom);

GetOptions("oid=s"=>\$org_id,
	   "v=s" => \$version,
	   "z|zoom=s" => \$zoom,
	   );


my $db = new CoGe::Genome;
my ($org) = $db->get_org_obj->search(organism_id=>$org_id);
my $chars = 10 * 2^$zoom;
foreach my $di ($org->data_information)
  {
    next if $version && $di->version ne $version;
    my $go = $di->genomic_sequences->next;
    my $chr_len = $db->get_genomic_sequence_obj->get_last_position($di);
    my $tot = ceil ($chr_len/$chars);
#    print "Total number of images to be generated: ", $chr_len/1280,"\n";
    my $count = 0;
    foreach (my $i=0; $i<= $chr_len; $i+=$chars)
      {
	$count++;
	my $cmd = "http://toxic.berkeley.edu/CoGe/tiler.pl?";
	$cmd .= "&x=$i";
	$cmd .= "&z=7";
	$cmd .= "&iw=256";
	$cmd .= "&di=".$di->id."";
	$cmd .= "&chr=".$go->chr."";
#	$cmd .= " o=$org_id ";
#	$cmd .= " v=$version ";
#	$cmd .= " | /dev/null";
	print "getting $cmd\n";
	print "\timage $count / $tot\n";
	get("$cmd");
      }
  }
