#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CoGeX;


my ($help, $dsgid, $version, $type, $chr, $DEBUG);

GetOptions ("h|help" =>  \$help,
	    "dsgid=i" => \$dsgid,
	    "t|type=s"    => \$type,
	    "d|debug"     => \$DEBUG,
            "c|chr=s"     => \$chr,
	    );

$help = 1 unless ($dsgid);
help() if $help;



my $coge = CoGeX->dbconnect();
my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);

print join ("\t", qw(TYPE CHR START STOP STRAND NAME)),"\n";

foreach my $ds ($dsg->datasets)
  {
    print STDERR "working on ",$ds->name,"\n";
    foreach my $feat ($ds->features)
      {
	next if $type && $feat->type->name ne $type;
	print join ("\t", $feat->type->name, $feat->chr, $feat->start, $feat->stop, $feat->strand, join (", ", $feat->names)),"\n";
      }
  }
		


sub help
  {
    print qq{
Welcome to $0

This program generates a table of all the features in a genome

Options:
  
  -h | -help               prints this message

  -dsgid                   CoGe dataset group database id

  -t | -type               specify the type of feature from which to extract the sequence.
                           E.g. gene, CDS, mRNA, tRNA, etc.

  -d | -debug              print debugging messages


};
exit(0);
  }
