#!/usr/bin/perl

use strict;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my %orgs;

foreach my $org ($coge->resultset('Organism')->all())
  {
    #let's take care of some name correction business
    if ($org->name =~ /^\s/ || $org->name =~ /\s$/)
      {
	my $name = $org->name;
	$name =~ s/^\s+//;
	$name =~ s/\s+$//;
	$org->name($name);
	$org->update;
      }
    if ($org->description =~ /^\s/ || $org->description =~ /\s$/)
      {
	my $description = $org->description;
	$description =~ s/^\s+//;
	$description =~ s/\s+$//;
	$org->description($description);
	$org->update;
      }
    if ($org->description =~ /\.$/)
      {
	my $description = $org->description;
	$description =~ s/\.$//;
	$org->description($description);
	$org->update;
      }

    if ($orgs{$org->name})
      {
	print $org->name,"\n";
	print $org->description,"\n";
	print $org->id,"\n";
	print "--------\n";
	print $orgs{$org->name}->name,"\n";
	print $orgs{$org->name}->description,"\n";
	print $orgs{$org->name}->id,"\n";
	print "\n\n";

	foreach my $dsg ($org->dataset_groups)
	  {
	    $dsg->organism_id($orgs{$org->name}->id);
	    $dsg->update;
	  }
	$org->delete;
      }
    else
      {
	$orgs{$org->name}=$org;
      }
  }
