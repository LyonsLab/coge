#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

use vars qw($coge $name $desc $help);

$coge = CoGeX->dbconnect();
GetOptions ("name|n=s"=>\$name,
	    "desc|d=s"=>\$desc,
	    "help|h"=>\$help,
	   );

if ($help)
  {
    help();
  }

my $search = {};
$search->{description}={like=>"%".$desc."%"} if $desc;
$search->{name}={like=>"%".$name."%"} if $name;

foreach my $org ($coge->resultset("Organism")->search($search))
  {
    print join ("\t", $org->name, $org->description),"\n";
  }

sub help
  {
    print qq{
Welcome to $0.
This program searches for organisms in CoGe matching a name and/or description
and generates a text list of the organism names and descriptions which match.
If nothing is specified, all organisms are listed.

Usage $0 -name <org_name> -desc <org_description>

Options:

 -name | -n    name of organisms to search (wildcards pre and post pended)

 -desc | -d    description of organisms to search (wildcards pre and post pended)

 -hehlp | -h   print this message

};
    exit;
  }
