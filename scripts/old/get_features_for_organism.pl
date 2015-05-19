#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CoGeX;

my ($help, $org, $version, $type, $chr, $DEBUG);

GetOptions ("h|help" =>  \$help,
	    "o|org=s" => \$org,
	    "v|version=s" => \$version,
	    "t|type=s"    => \$type,
	    "d|debug"     => \$DEBUG,
            "c|chr=s"     => \$chr,
	    );

$help = 1 unless ($org);
help() if $help;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

foreach my $ds ($coge->resultset('Dataset')->search({organism_id=>$org, version=>$version}))
  {
    print STDERR "working on ",$ds->name,"\n";
    foreach my $feat ($ds->features)
      {
	next if $type && $feat->type->name ne $type;
	print ">",join ("\t", sort $feat->names, "type: ".$feat->type->name),"\n";
	print $feat->genomic_sequence,"\n";
      }
  }

sub help
  {
    print qq{
Welcome to $0

This program generates a list of names for all features in an organism

Options:

  -h | -help               prints this message

  -o | -organism           specify database id for an organism

  -v | -version            specify the version of the dataset to use. If not specified
                           this will use the most current version of an organism

  -t | -type               specify the type of feature from which to extract the sequence.
                           E.g. gene, CDS, mRNA, tRNA, etc.

  -c | -chr                specify a chromosome of an organism.  IF not specified, all chrs
                           will be used.

  -d | -debug              print debugging messages

};
exit(0);
  }
