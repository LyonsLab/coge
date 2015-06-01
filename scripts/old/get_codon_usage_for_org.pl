#!/usr/bin/perl -w

use strict;
use CoGeX;

my $org = shift;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

($org) = $coge->resultset('Organism')->resolve($org);

print STDERR $org->name,"\n";
my ($code, $code_type);
my $header = 0;
foreach my $ds ($org->current_datasets)
  {
    foreach my $feat ($ds->features({"feature_type.name"=>"CDS"},{join=>"feature_type"}))
      {
	($code, $code_type) = $feat->genetic_code() unless $code;
	my %aa;
	grep {$aa{$_}=0} values %$code;
	unless ($header)
	  {
	    $header=1;
	    print join ("\t", (map {$_."%"} sort keys %$code),(map{$_."%"} sort keys %aa), "AT%", "GC%"),"\n";
	  }
	my $codon_total=0;
	my ($codon) = $feat->codon_frequency(counts=>1);
	grep {$codon_total+=$_} values %$codon;
	my $aa_total=0;
	foreach my $tri (keys %$code)
	  {
	    $aa{$code->{$tri}}+=$codon->{$tri};
	    $aa_total+=$codon->{$tri};
	  }
	my ($gc, $at) = $feat->gc_content();
	print join ("\t", (map {sprintf("%.2f", 100*$codon->{$_}/$codon_total)} sort keys %$code), (map {sprintf("%.2f",100*$aa{$_}/$aa_total) } sort keys %aa), 100*$at, 100*$gc ),"\n";
      }
  }
