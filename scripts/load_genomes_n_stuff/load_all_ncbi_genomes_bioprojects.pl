#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use Data::Dumper;
use XML::Simple;
use Parallel::ForkManager;

my $bpids = get_NCBI_bioprj();

my $pm = new Parallel::ForkManager(20);
my $out = "Load_all_NCBI.".sprintf( "%04d-%02d-%02d-%02d:%02d:%02d",
                 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime)).".log";

foreach my $bpid (@$bpids)
  {
    my $pid = $pm->start and next;
    my $output;
    $output.= "BioProject $bpid Started\n";
    if (check_bioprj_genome($bpid))
      {
	my $ntids = get_bioprj_nt($bpid);
	my @accns;
	foreach my $ntid (@$ntids)
	  {
	    my $accn = get_nt_accn($ntid);
	    push @accns, $accn;
	  }
	if (@accns)
	  {
	    my $processed_accns = process_accns(\@accns, \$output);
	    run_load_genome($processed_accns, \$output);
	  }
      }
    $output .= "BioProject $bpid Finished\n";
    open (OUT, ">>$out");
    print OUT $output,"\n\n";;
    close OUT;
    $pm->finish;
  }
$pm->wait_all_children;

sub get_NCBI_bioprj
  {
    my $esearchgenome = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=all%5Bfilter%5D&retmax=999999"; #get all genome ids
    my $entry = get($esearchgenome);
    my @bpids;
    while ($entry =~ /<Id>(\d+)<\/Id>/ig)
      {
	my $bpid = $1;
	push @bpids, $bpid;
      }
    return \@bpids;
  }

sub check_bioprj_genome
  {
    my $bpid = shift;
    my $esummary = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&retmode=text&complexity=0&id=";
    my $summary = get($esummary.$bpid);

    return 1 if $summary=~/Genome sequencing/;
    return 0;
  }

sub get_bioprj_nt
  {
    my $bpid = shift;
    my $elink = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=nucleotide&dbfrom=bioproject&id=";
    my $entry = get($elink.$bpid);
    my @ids;
    return \@ids unless $entry;
    my $xml = XMLin($entry);
    my $items = ref ($xml->{LinkSet}{LinkSetDb}{Link}) =~ /array/i ? [@{$xml->{LinkSet}{LinkSetDb}{Link}}] : [$xml->{LinkSet}{LinkSetDb}{Link}];
     foreach my $item (@$items)
       {
	 my $id = $item->{Id};
	 push @ids, $id if $id;
       }
     return \@ids;
  }

sub get_nt_accn
    {
      my $ntid = shift;
      my $esummary = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&retmode=text&complexity=0&id=";
      my $entry = get($esummary.$ntid);
      my ($accn) = $entry =~ /<Item Name="Caption" Type="String">(.+?)<\/Item>/;
      return $accn;
    }

sub process_accns
      {
	my $accns = shift;
	my $output = shift;
	$$output .= join ("\t", @$accns)."!\n";
	my $found=0;
	foreach my $accn (@$accns)
	  {
	    $found = $accn if $accn =~ /[A-Z]{4}0{8}/;
	  }
	if ($found)
	  {

	    $$output .= ("#"x20)."\n";
	    $$output .= $found."\n";
	    $$output .= ("#"x20)."\n";
	    return [$found];
	  }
	return $accns;
      }

sub run_load_genome
  {
    my $accns = shift;
    my $output = shift;
    my $prog = '/home/elyons/projects/CoGeX/scripts/load_genomes_n_stuff/genbank_genome_loader.pl';
#    $prog .= " -autoupdate";
    $prog .= " -td '/tmp/gb/'";
    $prog .= " -go";
    $prog .= " -autoskip";
    $prog .= " -delete_src_file";
    $prog .= " -max_entries 10000";
    my $run = $prog." -accn ".join (" -accn ", @$accns);
    $$output .=  $run."\n";
    open (IN, $run." |");
    while (<IN>)
      {
#	$previously_loaded =1 if /previously loaded/;
#	$skipped=1 if /skipping/;
	$$output .= $_;
      }
    close IN;
  }
