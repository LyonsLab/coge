#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
my $infile;
my $outfile1;
my $outfile2;

GetOptions (
	    "infile=s" => \$infile,
	    "outfile1=s"    => \$outfile1,
	    "outfile2=s"    => \$outfile2,
	   );

my $input =[];
open (IN, $infile)  || die "Can't open $infile for reading: $!";
while (<IN>)
  {
    chomp;
    next unless $_;
    push @$input, $_;
  }
close IN;

$input = convert_blast_genomic_names(data=>$input, outfile=>$infile.".new") if ($infile =~ /genomic/);

my %seen1 = ();
my %seen2 = ();
my @output1;
my @output2;

foreach my $item (@$input)
  {
    my @line = split/\t/, $item;
    unless ($line[0] && $line[1])
      {
       print STDERR "Skipping line because it is undefined: $_\n";
       next;
      }
    unless ($line[0]=~/\|\|/ && $line[1]=~/\|\|/)
      {
	print STDERR "Skipping entry because it is not formatted correctly with '||':\nEntry 1: $line[0]\nEntry 2: $line[1]\n";
	next;
      }

    my @item1 = split/\|\|/, $line[0];
    my @item2 = split/\|\|/, $line[1];
    unless (defined $item1[0] && defined $item1[1] && defined $item1[2])
      {
	print STDERR "Skipping printing output because of missing values:\n";
	print STDERR join ("\n", $item1[0], $item1[1], $item1[2]),"\n";
	next;
      }
    unless (defined $item2[0] && defined $item2[1] && defined $item2[2])
      {
	print STDERR "Skipping printing output because of missing values:\n";
	print STDERR join ("\n", $item2[0], $item2[1], $item2[2]),"\n";
	next;
      }
    push @output1, join ("\t", $item1[0], $item1[1], $item1[2], $line[0]) unless $seen1{$line[0]};
    push @output2, join ("\t", $item2[0], $item2[1], $item2[2], $line[1]) unless $seen2{$line[1]};
    $seen1{$line[0]} =1;
    $seen2{$line[1]} =1;
  }
open (OUT1, ">$outfile1") || die "Can't open $outfile1 for writing: $!";
print OUT1 join ("\n", @output1);
close OUT1;

open (OUT2, ">$outfile2") || die "Can't open $outfile2 for writing: $!";
print OUT2 join ("\n", @output2);
close OUT2;

sub convert_blast_genomic_names
  {
    my %opts = @_;
    my $data = $opts{data};
    my $outfile = $opts{outfile};

    my ($hits_order1, $hits_order2) = order_hits($data);

    my @output;
    foreach my $item (@$data)
      {
	my @line = split/\t/, $item;
	my @item1 = split/\|\|/, $line[0];
	my @item2 = split/\|\|/, $line[1];
	#add start and stop positions to items if not present (e.g. genomci sequence hits)
	my ($ori1, $ori2); #strand/orientation of features/hits
	$ori1 = $item1[4];
	$ori2 = $item2[4];
	unless (defined $item1[1] && defined $item1[2] && defined $item1[4])
	  {
	    $item1[1] = $line[6] unless defined $item1[1];
	    $item1[2] = $line[7] unless defined $item1[2];
	    ($item1[1], $item1[2]) = ($item1[2], $item1[1]) if $item1[1]>$item1[2];
	    #determine orientation of the hit relative to the gene (if present)
	    $ori1 = $line[6] > $line[7] ? -1 : 1;
	  }
	unless (defined $item2[1] && defined $item2[2] && defined $item2[4])
	  {
	    $item2[1] = $line[8] unless defined $item2[1];
	    $item2[2] = $line[9] unless defined $item2[2];
	    ($item2[1], $item2[2]) = ($item2[2], $item2[1]) if $item2[1]>$item2[2];
	    $ori2 = $line[8] > $line[9] ? -1 : 1;
	  }
	my $ori = $ori1 eq $ori2 ? 1 : -1;
	#for genomic sequences -- needto add genomic_hit and remove 'gi|' or 'lcl|' from begining of fasta header name
	unless ($item1[3])
	  {
	    $item1[3] = "genomic_hit";
	    $item1[0] =~ s/^gi\|//;
	    $item1[0] =~ s/^lcl\|//;
	  }
	unless ($item2[3])
	  {
	    $item2[3] = "genomic_hit";
	    $item2[0] =~ s/^gi\|//;
	    $item2[0] =~ s/^lcl\|//;
	  }
	unless ($item1[4])
	  {
	    $ori1 = $ori eq 1 ? $ori2 : $ori2*-1;
	  }
	unless ($item2[4])
	  {
	    $ori2 = $ori eq 1 ? $ori1 : $ori1*-1;
	  }
	$item1[4] = $ori1 unless $item1[4];
	$item2[4] = $ori2 unless $item2[4];
	unless ($item1[5])
	  {
	    $item1[5] = "genomic_hit";
	    $item1[6] = "0";
	    my $key = $item1[0]."-".$line[6]."-".$line[7];
	    $item1[7] = $hits_order1->{$key};
	  }

	unless ($item2[5])
	  {
	    $item2[5] = "genomic_hit";
	    $item2[6] = "0";
	    my $key = $item2[0]."-".$line[8]."-".$line[9];
	    $item2[7] = $hits_order2->{$key};
	  }
	my $name1 = join ("||", @item1);
	my $name2 = join ("||", @item2);
	$line[0] = $name1;
	$line[1] = $name2;
	push @output, join ("\t", @line);
      }
    open (OUT, ">$outfile") || die "Can't open $outfile for writing: $!";
    print OUT join "\n", @output;
    close OUT;
    return \@output;
  }

sub order_hits
  {
    my $hits = shift;
    my (%data1, %data2); #data to be returned
    my (@to_sort1, @to_sort2); #data to be sorted
    foreach my $item (@$hits)
      {
	my @items = split/\t/, $item;
	#'||' is the character used to delimit genomic feature information, if not present, it is a genomic hit which needs sorting
	unless ($items[0] =~ /\/\//)
	  {
	    $items[0] =~ s/gi\|//;
	    $items[0] =~ s/^lcl\|//;
	    push @to_sort1, [$items[0],$items[6], $items[7]];
	  }
	unless ($items[2] =~ /\/\//)
	  {
	    $items[1] =~ s/gi\|//;
	    $items[1] =~ s/^lcl\|//;
	    push @to_sort2, [$items[1],$items[8], $items[9]];
	  }
      }
    my $sorted1 = sort_hits(\@to_sort1) if @to_sort1;
    my $sorted2 = sort_hits(\@to_sort2) if @to_sort2;
    return $sorted1, $sorted2;
  }

sub sort_hits
  {
    my $data = shift;
    my %sorted;
    my $count =1;
    foreach my $item (sort {
      $a->[0] cmp $b->[0] ||
      $a->[1] <=> $b->[1] ||
      $a->[2] <=> $b->[2]
			   } @$data)
      {
       my $key = $item->[0]."-".$item->[1]."-".$item->[2];
       next if $sorted{$key};
       $sorted{$key}=$count;
       $count++;
      }
    return \%sorted;
  }
