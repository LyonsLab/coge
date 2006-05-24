#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGe::Genome;
use LWP::Simple;

my ($org_id, $version);

GetOptions("oid=s"=>\$org_id,
	   "v=s" => \$version,
	   );

my $db = new CoGe::Genome;
my ($org) = $db->get_org_obj->search(organism_id=>$org_id);
foreach my $di ($org->data_information)
  {
    next if $version && $di->version ne $version;
    my $go = $di->genomic_sequences->next;
    my $chr_len = $db->get_genomic_sequence_obj->get_last_position($di);
    
    foreach (my $i=0; $i<= $chr_len; $i+=1280)
      {
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
	get("$cmd");
      }
  }
