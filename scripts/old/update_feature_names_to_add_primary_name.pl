#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my ($help, $debug, $orgname, $name_string, $DEBUG, $go);

GetOptions ("h|help" =>  \$help,
            "o|org=s" => \$orgname,
            "d|debug"     => \$DEBUG,
            "n|name=s" =>\$name_string,
	    "go=s"=>\$go,
            );
$| =1;
$go = 1 unless defined $go;
unless ($orgname && $name_string)
  {
    print qq{
Usage: $0 -o <org name or id> -n <regular expression name search>

This program sets feature names to primary_names if the name matches the regular expression.  E.g.:
};
print "\t\$ ",$0,' -o Arabidopsis -n "^at\dg\d\d\d\d\d$"';
print qq {
which will match any name like at1g01010
};
    exit;
  }

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my ($org) = $coge->resultset('Organism')->resolve($orgname);

print "working on ".$org->name.": ".$org->description."\n";

foreach my $ds ($org->datasets)
  {
    print "\tdataset: ",$ds->name,": ",$ds->description, " (v".$ds->version.")"."\n\tProgress (each '.' is 100 names updated):";
    my $count = 0;
    foreach my $feat ($ds->features)
      {
	foreach my $name ($feat->feature_names)
	  {
	    if ($name->name=~/$name_string/i)
	      {
		print $name->name, " matches ", $name_string,"\n" if $DEBUG;
		print "." unless $count %100;
		$count++;
		$name->update({primary_name=>1}) if $go;
	      }
	  }
      }
    print "\n";
  }
