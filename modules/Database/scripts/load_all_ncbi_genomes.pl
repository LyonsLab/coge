#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Long;
use XML::Simple;

$| =1;

my ($DEBUG, $forks, $pause, $autoupdate);

GetOptions (
	    "debug"=>\$DEBUG,
	    "forks|f=i"=>\$forks,
	    "pause|p"=>\$pause,
	    "autoupdate"=>\$autoupdate,
	    );
$forks = 1 unless defined $forks;


#print Dumper get_genomes_for_genomeprj(116);
#exit;


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
# 9 == chromosomes (complete and WGS) and plasmids (NC_ NS_ NZ_)

my $base_url = "http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?";
#taxid type
my @taxids =(
#	     [2157, 9], #archaea, get WGS
	     [2,    9], #bacteria, get WGS
#	     [2759, 0], #euks, not ready for WGS -- lots of data with minimal annotations.  Will get there though!
#	     [10239,9], #viruses, phages, 
#	     [12884,9], #viroids
	    );

my $pm = new Parallel::ForkManager($forks);
# foreach my $item (@taxids)
#   {
#     my ($id, $type) = @$item;
#     my $url = $base_url."type=$type"."&taxid=$id";
#     print "Fetching $url\n";
#     my $content =get($url);
#     my @tables = split /<\/table>/, $content;
#     my @rows = split /<tr>/,$tables[8];
#     my %data;
#     foreach my $row (@rows)
#       {

# 	next unless $row =~/^<td>/;
# 	$row =~ s/\n|\r//g;
# 	my @cols = split /<td.*?>/, $row;
# 	foreach my $item (@cols)
# 	  {
# 	    $item =~ s/<.*?>//g;
# 	  }
# 	next unless $cols[3];
# 	next if $cols[3] =~ /accession/;
# 	$cols[2] =~ s/&nbsp;//g;
# 	$cols[2] =~ s/chromosome//i;
# 	$cols[2] =~ s/chr//i;
# 	$cols[2] =~ s/^\s+//;
# 	$cols[2] =~ s/\s+$//;
# 	$cols[2] = undef if $cols[2] && ($cols[2] eq " " || $cols[2] =~ /nbsp/);
# 	$cols[2] = "unknown" unless $cols[2];
# 	#$data{organism}{chromosome} = [accns]
# 	push @{$data{$cols[1]}{$cols[2]}}, $cols[3];

#       }

#     my $orgs = process_orgs(\%data);

my $genomes = get_NCBI_genomes();


my $prog = '/home/elyons/projects/CoGeX/scripts/genbank_genome_loader.pl';
$prog .= " -autoupdate" if $autoupdate;

foreach my $item (@$genomes)
  {
    sleep 1;
    #need to break apart organelle genomes from nuclear genomes
    my @tmp; #temp place to store genomes in case some need to be processed separately in addition to the whole genome.  Such as is the case with mitochondria, chloroplasts, plastids.
    if (@{$item->{orgs}} > 1) #don't need to check if there is only one item for the genome
      {
	my @subgenome_names = qw(mitochondrion mitochondrial chloroplast apicoplast plastid);
	my $i=0; #array index;
	my %no_organelles= (org=>[], accns=>[], gids=>[]); #storage for genomes without organelles
	foreach my $org (@{$item->{orgs}})
	  {
	    my $organelle_found=0;
	    foreach my $check (@subgenome_names)
	      {
		if ($org =~ /$check/)
		  {
		    $organelle_found=1;
		    #add organelle as separate genome
		    push @tmp, {
				org=>[$org],
				accns=>[$item->{accns}[$i]],
				gids=>[$item->{gids}[$i]],
			       };
		  }
	      }
	    unless ($organelle_found)
	      {
		push @{$no_organelles{org}}, $org;
		push @{$no_organelles{accns}}, $item->{accns}[$i];
		push @{$no_organelles{gids}}, $item->{gids}[$i];
	      }
	    $i++;
	  }
	push @tmp, \%no_organelles;
      }
    else #only one genome
      {
	push @tmp, $item;
      }

    foreach my $genome (@tmp)
      {
	#	$pm->start and next;
	print Dumper $genome;
	my $run = $prog." -accn ".join (" -accn ", @{$genome->{accns}});
	#$run .= " -chr '$cols[2]'" if $cols[2];
	$run .= " -td '/tmp/gb/'";
	$run .= " -go";
	print "\n",$run,"\n";
	my $previously_loaded = 0;
#	next;
	open (IN, $run." |");
	while (<IN>)
	  {
	    $previously_loaded =1 if /previously loaded/;
	    print $_;
	  }
	close IN;
	if ($pause && !$previously_loaded)
	  {
	    print "Paused:  press enter to continue\n";
	    <STDIN>;
	  }
	#	$pm->finish;
      }
  }
$pm->wait_all_children;


sub process_orgs
  {
    my $data = shift;
    my %orgs;
    foreach my $org(keys %$data)
      {
	foreach my $chr (keys %{$data->{$org}})
	  {
	    if (scalar @{$data->{$org}{$chr}} > 1 )
	      {
		print Dumper $data->{$org}{$chr},"\n";
		foreach my $accn (@{$data->{$org}{$chr}})
		  {
		    my $project_id = get_project($accn);
		    push @{$orgs{$org}{$project_id}}, {accn=>$accn};
		  }
	      }
	    else
	      {
		push @{$orgs{$org}{"NA"}}, {accn=>$data->{$org}{$chr}[0]};
	      }
	  }
      }
    return \%orgs;
  }

sub get_project
  {
    my $accn = shift;
    $accn = join ",", @$accn if ref ($accn) =~ /array/i;
    my $summary_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gbwithparts&retmode=text&complexity=3&id=";
    my $url = $summary_url.$accn;
#    print $url,"\n";
    my $summary = get($url);
    unless ($summary)
      {
	$summary = get($url);
      }
    unless ($summary)
      {
	print "Error getting $url\n";
	next;
      }
	  
    my %data;
    while ($summary =~/ACCESSION\s+(\S+).*?DBLINK\s+Project:(\d+)/gxs)
      {
	$data{$1} = $2;

      }
#    print scalar keys %data,"\n";
    return \%data;
  }


sub get_NCBI_genomes_old
  {
    my $url = "ftp://ftp.ncbi.nih.gov/genomes/genomeprj/gp.xml";
    my $file = "/home/elyons/tmp/gb_genomes.xml";
    getstore ($url, $file) unless -r $file && -A $file > 2; #get the file if older than 2 days
    my $gpxml;
    open (IN, $file);
    while (<IN>)
      {
	$gpxml .= $_;
      }
    $gpxml =~ s/<\?xml.*?>//;
    $gpxml =~ s/<!DOCTYPE.*?>//;
    $gpxml=~s/gp://g;
    my %genomes;
    my $count = 0;
    while ($gpxml =~ /ProjectID>(\d+)<\/ProjectID>/g)
      {
	my $gid = $1;
#	print $gid,"\n";
	my ($accns, $orgs) = get_genome_accns($gid);
	$genomes{$gid} = {accns=>$accns, orgs=>$orgs};# if @$accns;
	$count++;
#	last if $count == 100;
      }
    return \%genomes;
  }

sub get_NCBI_genomes
  {
    my $esearchgenomeprj = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genomeprj&term=all%5Bfilter%5D&retmax=999999"; #genomeprj let's me assemble genomes when there are multiple assemblies (no viruses)
    my $esearchgenome = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genomeprj&term=all%5Bfilter%5D&retmax=999999"; #get all genome ids
    my @genomes;


    #get genome projects
    my $entry = get($esearchgenomeprj);
    my %projects; #store genome_id=genomeprj_id
    my $i =0;
    while ($entry =~ /<Id>(\d+)<\/Id>/ig)
      {
#	last if $i == 2000;
	my $gpid = $1;
	my $gids = get_genomes_for_genomeprj($gpid); #metagenomes, etc won't have genomes. . . 
	my @gids;
	foreach my $gid (@$gids)
	  {
	    next unless $gid; #many genome projects have no specific genomes -- metagenomes, ENCODE
	    $projects{$gid}=$gpid;
	    push @gids, $gid;
	  }
	push @genomes, {gids=>\@gids,
		       gpid=>$gpid} if @gids;
	$i++;
      }
    #get genomes
    $entry = get($esearchgenome);
    while ($entry =~ /<Id>(\d+)<\/Id>/ig)
      {
	next;
	my $gid = $1;
	next if $projects{$gid};
	push @genomes, {gids=>[$gid]};
      }
    $i =0;
    foreach my $genome (@genomes)
      {
#	last if $i ==50;
	my $gids = join ",", @{$genome->{gids}};
	my ($accns, $orgs) = get_genome_accns($gids);
	$genome->{accns}=$accns;
	$genome->{orgs}=$orgs;# if @$accns;
	$i++;
      }
    return \@genomes;
  }

sub get_genomes_for_genomeprj
    {
      my $gid = shift;
      my $elink = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=genome&dbfrom=genomeprj&id=";
      my $entry = get($elink."$gid");
      my $xml = XMLin($entry);
      my $items = ref ($xml->{LinkSet}{LinkSetDb}{Link}) =~ /array/i ? [@{$xml->{LinkSet}{LinkSetDb}{Link}}] : [$xml->{LinkSet}{LinkSetDb}{Link}];
      my @ids;
      foreach my $item (@$items)
	{
	  my $id = $item->{Id};
	  push @ids, $id if $id;
	}
      return \@ids;
    }

sub get_genome_accns
    {
      my $gid = shift;
      my $esumg = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=genome&retmode=text&complexity=0&id=";
      my $item = get($esumg.$gid);
      my @accns;
      my @orgs;
      while ($item=~ /<Item Name="Caption" Type="String">(.*?)<\/Item>/g)
	{
	  push @accns, $1;
	}
      while ($item=~ /<Item Name="Title" Type="String">(.*?)<\/Item>/g)
	{
	  print STDERR $1,"\n";
	  push @orgs, $1;
	}
      return \@accns, \@orgs;
    }

sub get_genome_project_id_for_accn
   {
     my $accn = shift;
     my $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=";
     my $elink = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=genomeprj&id=";
     my $entry = get ($esearch.$accn);
     my ($id) = $entry =~ /<Id>(\d+)<\/Id>/;
     my $entry2 = get ($elink.$id);
     my $gid;
     while ($entry2 =~/<Id>(\d+)<\/Id>/)
       {
	 $gid = $1 if $1 ne $id;
       }
     return $gid;
   }
