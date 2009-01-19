#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use Data::Dumper;

my ($DEBUG, $oid1, $oid2, $input_file);

GetOptions(
	   "DEBUG"=>\$DEBUG,
	   "oid1|org1|o1=s"=>\$oid1,
	   "oid2|org2|o2=s"=>\$oid2,
	   "input|file|i|f=s"=>\$input_file,
	  );


my $coge = CoGeX->dbconnect();
my $org1 = $coge->resultset('Organism')->resolve($oid1);
my $org2 = $coge->resultset('Organism')->resolve($oid2);
my $order1 = get_gene_order($org1);
my $order2 = get_gene_order($org2);
#print Dumper [keys %$order2];
#print map {$_."\n"} sort {$order2->{'contig_NZ_AAKB02000001'}{$a} <=> $order2->{'contig_NZ_AAKB02000001'}{$b}} keys %{$order2->{'contig_NZ_AAKB02000001'}};
#print Dumper $order2;
open (IN, $input_file);
while (<IN>)
  {
    chomp;
    my @line = split/\t/;
    my @item1 = split/\|\|/, $line[1];
    my @item2 = split/\|\|/, $line[5];
    $line[2] = $order1->{$item1[0]}{$item1[3]};
    $line[3] = $order1->{$item1[0]}{$item1[3]};
    $line[6] = $order2->{$item2[0]}{$item2[3]};
    $line[7] = $order2->{$item2[0]}{$item2[3]};
    print join "\t", @line,"\n";
  }
close IN;




sub get_gene_order
  {
    my $org = shift;
    my %data;
    foreach my $ds ($org->current_datasets())
      {
	my %chr;
	foreach my $feat (sort {$a->start <=> $b->start} $coge->resultset('Feature')->search({
											      datAset_id=>$ds->id,
											      feature_type_id=>3}))
	  {
	    $chr{$feat->chromosome}++;
	    my ($name) = $feat->names;
	    $data{$feat->chromosome}{$name} = $chr{$feat->chromosome};
	  }
      }
    return \%data;
  }





