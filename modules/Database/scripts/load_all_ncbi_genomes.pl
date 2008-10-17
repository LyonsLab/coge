#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use Data::Dumper;
use Parallel::ForkManager;

# for understanding accession nomeclature from NCBI see http://www.ncbi.nlm.nih.gov/RefSeq/key.html#accessions
# genbank genlist.cgi type codes:
# 0 == Chromosomes and plasmids (NC_)
# 1 == Chromosomes (NC_)
# 2 == plasmids (NC_)
# 3 == Whole shotgun sequences (NZ_)
# 4 == mito, plastids, plasmids, nucleomorphs
# 5 == viruses
# 6 == phages
# 7 == all plasmids for all proks? (NC_)
# 8 == Chromosomes (NC_ and NS_)
# 9 == chromosomes and plasmids (NC_ NS_ NZ_)
my @pages = qw(
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2157&type=1&name=Archaea%20Complete%20Chromosomes
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2&type=0
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2759&type=8&name=Eukaryotae%20Complete%20Chromosomes
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=10239&type=5&name=Viruses
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=10239&type=6&name=Phages
http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=12884&type=0&name=Viroids
);
my $base_url = "http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?";

#taxids
#2157 == archaea
#2    == bacteria
#2759 == euks

#taxid type
my @taxids =(
	     [2157, 9], #archaea, get WGS
	     [2,    9], #bacteria, get WGS
	     [2759, 0], #euks, not ready for WGS -- lots of data with minimal annotations.  Will get there though!
	     [10239,9], #viruses, phages, 
	     [12884,9], #viroids
	    );

my $pm = new Parallel::ForkManager(20);
foreach my $item (@taxids)
  {
    my ($id, $type) = @$item;
    my $url = $base_url."type=$type"."&taxid=$id";
    my $content =get($url);
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
	$run .= " -td '/tmp/gb/'";
	$run .= " -go";
	print "\n",$run,"\n";
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
