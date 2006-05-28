#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGe::Genome;
use LWP::Simple;
use POSIX;
use Benchmark;

my ($org_id, $version, $url);

GetOptions("oid=s"=>\$org_id,
	   "v=s" => \$version,
#	   "z|zoom=s" => \$zoom,
	   "u|url=s" => \$url
	   );


$url = "synteny.cnr.berkeley.edu" unless $url;
$url = "http://" .$url unless $url =~ /http:\/\//;
$url .= "/CoGe/tiler.pl";
my $db = new CoGe::Genome;
my ($org) = $db->get_org_obj->search(organism_id=>$org_id);

foreach my $di ($org->data_information)
  {
    next if $version && $di->version ne $version;
    my $go = $di->genomic_sequences->next;
    next unless $go;
    my $chr_len = $go->get_last_position($di);
    next unless $chr_len;
    my $max_zoom = ceil (log10($chr_len)/(log10(10)*log10(2)));
    foreach my $zoom (0..$max_zoom)
      {
	my $chars = 10 * 2**$zoom;
	my $tot = ceil ($chr_len/$chars);
	print "Total number of bp: $chr_len\n";
	print "$chars characters per tile at zoom level $zoom\n";
	print "Total number of images to be generated: ", $tot,"\n";
#	next;
	my $count = 0;
	foreach (my $i=0; $i<= $chr_len; $i+=$chars)
	  {
	    $count++;
	    my $cmd = "$url"."?";
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
  }
