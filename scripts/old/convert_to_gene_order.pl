#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use Getopt::Long;
use Data::Dumper;

my ($DEBUG, $dsgid1, $dsgid2, $input_file, $ftid1, $ftid2, $conffile);

GetOptions(
	   "DEBUG"=>\$DEBUG,
	   "dsgid1|dsg1|d1=s"=>\$dsgid1,
	   "dsgid2|dsg2|d2=s"=>\$dsgid2,
	   "input|file|i|f=s"=>\$input_file,
	   "ftid1|ft1=i"=>\$ftid1,
	   "ftid2|ft2=i"=>\$ftid2,
	   "config_file|cf=s"=>\$conffile,
	  );

my $P = CoGe::Accessory::Web::get_defaults($conffile);
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

my $order1 = get_gene_order(dsgid=>$dsgid1, ftid=>$ftid1) if $ftid1 == 3;
my $order2 = $dsgid1 == $dsgid2 && $ftid1 == $ftid2 ? $order1 : get_gene_order(dsgid=>$dsgid2, ftid=>$ftid2) if $ftid2 == 3;

$order1 = get_order_from_input(file=>$input_file, set=>1) unless $order1;
$order2 = get_order_from_input(file=>$input_file, set=>2) unless $order2;
#print Dumper $order1, $order2;
my %seen;
open (IN, $input_file);
while (<IN>)
  {
    chomp;
    my @line = split/\t/;
    my @item1 = split/\|\|/, $line[1];
    my @item2 = split/\|\|/, $line[5];
#    print $item1[0],"\n",$item1[3],"\n";
    if ($item1[4])
      {
	next unless $order1->{$item1[0]}{$item1[3]};
	$line[2] = $order1->{$item1[0]}{$item1[3]};
	$line[3] = $order1->{$item1[0]}{$item1[3]};
      }
    else
      {
	$line[2] = $order1->{$item1[0]}{$item1[1]}{$item1[2]};
	$line[2] = $order1->{$item1[0]}{$item1[2]}{$item1[1]} unless $line[2];
	$line[3] = $order1->{$item1[0]}{$item1[1]}{$item1[2]};
	$line[3] = $order1->{$item1[0]}{$item1[2]}{$item1[1]} unless $line[3];
      }
    if ($item2[4])
      {
	next unless $order2->{$item2[0]}{$item2[3]};
	$line[6] = $order2->{$item2[0]}{$item2[3]};
	$line[7] = $order2->{$item2[0]}{$item2[3]};
      }
    else
      {
	$line[6] = $order2->{$item2[0]}{$item2[1]}{$item2[2]};
	$line[6] = $order2->{$item2[0]}{$item2[2]}{$item2[1]} unless $line[6];
	$line[7] = $order2->{$item2[0]}{$item2[1]}{$item2[2]};
	$line[7] = $order2->{$item2[0]}{$item2[2]}{$item2[1]} unless $line[7];
      }
    print join "\t", @line,"\n" unless $seen{$line[2]}{$line[6]};
    $seen{$line[2]}{$line[6]}=1;
  }
close IN;

sub get_order_from_input
  {
    my %opts = @_;
    my $file = $opts{file};
    my $set = $opts{set};
    my %data;
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
	my @line = split/\t/;
	my @item = $set == 1 ? split(/\|\|/,$line[1]) : split/\|\|/, $line[5];
	my ($start, $stop) = sort {$a <=> $b} ($item[1],$item[2]);
	$data{$item[0]}{$start}{$stop}=1;
      }
    close IN;
    foreach my $chr (keys %data)
      {
	my $count =1;
	foreach my $start (sort {$a <=> $b} keys %{$data{$chr}})
	  {
	    foreach my $stop (sort {$a <=> $b} keys %{$data{$chr}{$start}})
	      {
		$data{$chr}{$start}{$stop}=$count;
		$count++;
	      }
	  }
      }
    return \%data;
  }

sub get_gene_order
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $ftid = $opts{ftid};
    $ftid = 3 unless defined $ftid; #default for CDS
    my %data;
    my %chr;
    foreach my $feat (sort {$a->chromosome cmp $b->chromosome || $a->start <=> $b->start}
		      $coge->resultset('Feature')->search(
							  {
							   feature_type_id=>[3,4,7],
							   dataset_group_id=>$dsgid
							  },{
							     join=>[{dataset=>'dataset_connectors'},'feature_names'],
							     prefetch=>['feature_names']}
							 ))
      {
	$chr{$feat->chromosome}++;
	foreach my $name($feat->names)
	  {
	    $data{$feat->chromosome}{$name} = $chr{$feat->chromosome};
	  }
      }
    return 0 unless keys %data;
    return \%data;
  }
