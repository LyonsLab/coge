#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGe::Genome;
use LWP::Simple;
use POSIX;
use Benchmark;

my ($org_id, $version, $zoom);

GetOptions("oid=s"=>\$org_id,
	   "v=s" => \$version,
	   "z|zoom=s" => \$zoom,
	   );


my $db = new CoGe::Genome;
my ($org) = $db->get_org_obj->search(organism_id=>$org_id);
my $chars = 10 * 2**$zoom;
foreach my $di ($org->data_information)
  {
    next if $version && $di->version ne $version;
    my $go = $di->genomic_sequences->next;
    my $chr_len = $db->get_genomic_sequence_obj->get_last_position($di);
    my $tot = ceil ($chr_len/$chars);
#    print "Total number of bp: $chr_len\n";
#    print "$chars characters per tile at zoom level $zoom\n";
#    print "Total number of images to be generated: ", $tot,"\n";
#    exit;
    my $count = 0;
    foreach (my $i=0; $i<= $chr_len; $i+=$chars)
      {
	$count++;
	my $cmd = "http://toxic.berkeley.edu/CoGe/tiler.pl?";
	$cmd .= "&x=$i";
	$cmd .= "&z=$zoom";
	$cmd .= "&iw=256";
	$cmd .= "&di=".$di->id."";
	$cmd .= "&chr=".$go->chr."";
#	$cmd .= " o=$org_id ";
#	$cmd .= " v=$version ";
#	$cmd .= " | /dev/null";
	print "getting $cmd\n";
	print "\timage $count / $tot\n";
	my $t0 = new Benchmark;
	get("$cmd");
	my $t1 = new Benchmark;
	my $time = timestr(timediff($t1, $t0));
	print "\tImage generation took $time\n";
      }
  }
