#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my %data;
foreach my $ds (sort {$b->id <=> $a->id} $coge->resultset('Dataset')->all())
  {
    next unless $ds->name =~ /gbk/;
    if ($data{$ds->name}{$ds->version})
      {
#	$ds->delete() unless ($ds->organism);
#	next unless ($ds->organism->id ne $data{$ds->name}{$ds->version}->organism->id);
	print $ds->name. " already exists:\n";
	print "\t", $ds->description,": ".$ds->id,"\n";
#	print "\t\t", $ds->organism->name,": ",$ds->organism->id,"\n";
	print "\t", $data{$ds->name}{$ds->version}->description,": ".$data{$ds->name}{$ds->version}->id,,"\n";
#	print "\t\t", $data{$ds->name}{$ds->version}->organism->name,": ",$data{$ds->name}{$ds->version}->organism->id,"\n";
	print "Deleting ".$ds->id."\n";
#	$ds->delete();
	next;
      }
    $data{$ds->name}{$ds->version}=$ds;

  }
