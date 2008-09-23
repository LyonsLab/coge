#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use Data::Dumper;
use Parallel::ForkManager;

my @pages = qw(
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2157&type=1&name=Archaea%20Complete%20Chromosomes
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2&type=1&name=Bacteria%20Complete%20Chromosomes
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2759&type=8&name=Eukaryotae%20Complete%20Chromosomes
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=10239&type=5&name=Viruses
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=10239&type=6&name=Phages
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=12884&type=0&name=Viroids
);


my $pm = new Parallel::ForkManager(20);
foreach my $page (@pages)
  {
    my $content =get($page);
    my @tables = split /<\/table>/, $content;
    my @rows = split /<tr>/,$tables[8];
    my $prog = '/home/elyons/projects/CoGeX/scripts/genbank_genome_loader.pl';
    foreach my $row (@rows)
      {
	$pm->start and next;
	next unless $row =~/^<td>/;
	$row =~ s/\n|\r//g;
	my @cols = split /<td.*?>/, $row;
	#    print Dumper \@cols;
	foreach my $item (@cols)
	  {
	    $item =~ s/<.*?>//g;
	  }
	next unless $cols[3];
	next if $cols[3] =~ /accession/;
	my $run = $prog." -accn $cols[3]";
	$cols[2] =~ s/chromosome//i;
	$cols[2] =~ s/chr//i;
	$cols[2] =~ s/^\s+//;
	$cols[2] =~ s/\s+$//;
	$cols[2] = undef if $cols[2] && ($cols[2] eq " " || $cols[2] =~ /nbsp/);
	$run .= " -chr '$cols[2]'" if $cols[2];
	$run .= " -go -td '/tmp/gb/'";
	print $run,"\n";
	open (IN, $run." |");
	while (<IN>)
	  {
	    print $_;
	  }
	close IN;
	$pm->finish;
      }
  }
$pm->wait_all_children;
