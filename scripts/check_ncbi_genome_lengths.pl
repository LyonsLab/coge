#!/usr/bin/perl -w

use strict;
use CoGeX;
use LWP::Simple;
use Data::Dumper;
use CoGe::Accessory::GenBank;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

foreach my $ds ($coge->resultset('Dataset')->search({data_source_id=>21}))
  {
    my %chrs;
    my %lengths;
    foreach my $feat ($ds->features({feature_type_id=>301}))
      {
	$chrs{$feat->chromosome}=$feat->stop;
	$lengths{$feat->stop}=1;
#	print join ("\t", $feat->chromosome, $feat->start, $feat->stop),"\n";
      }
    my $total =0;
    map{$total+=$_} keys %lengths;
    $lengths{$total}=1;
    $chrs{total}=$total;
    my $link = $ds->link;
    $link .= "&sendto=t";
    $link =~ s/dopt=gbwithparts&//;
#    print $link,"\n";
    my $content = get($link);
#    print $content,"\n";
#    while ($content =~ /REFERENCE\s+(.*?)\n/g)
#      {
#	print $1,"\n";
#      }
    my ($source) = $content =~ /source\s+(.*?)\n/;
    my ($length) = $source =~ /\.\.(\d+)/;
    sleep(1);

    print "No length extracted for source: $source\n" unless $length;
    next if $lengths{$length};
    print $ds->name, "\t", $length,"\n";
    map{print $_."\t".$chrs{$_}."\n"} keys %chrs;
    print "\n";
  }
