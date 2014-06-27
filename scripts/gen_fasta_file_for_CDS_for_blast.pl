#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;

my $coge = CoGeX->dbconnect();
my $dsgid=shift;
unless ($dsgid)
  {
    print qq{
Usage: $0 <dsg_id>
};
  }

my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
foreach my $ds ($dsg->datasets)
  {
    foreach my $feat ($ds->features({feature_type_id=>3}))
      {
	my ($chr) = $feat->chromosome;#=~/(\d+)/;
	my $name;
	foreach my $n ($feat->names)
	  {
	    $name = $n;
	    last unless $name =~ /\s/;
	  }
	$name =~ s/\s+/_/g;
	my $title = join ("||",$chr, $feat->start, $feat->stop, $name, $feat->strand, "CDS", $feat->id);
	my $seq = $feat->genomic_sequence(dsgid=>$dsg);
	next unless $seq;
	#skip sequences that are only 'x' | 'n';
	next unless $seq =~ /[^x|n]/i;
	print ">".$title."\n";
	print $seq,"\n";
      }
  }
