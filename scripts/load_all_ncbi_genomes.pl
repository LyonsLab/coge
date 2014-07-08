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
my @skipped; #array to hold commands that had skipped reloading.  Will be printed out at end of program for user to manually go through;

#print Dumper get_genomes_for_genomeprj(116);
#exit;

# for understanding accession nomeclature from NCBI see:
# all: http://www.ncbi.nlm.nih.gov/Sequin/acc.html
# refseq: http://www.ncbi.nlm.nih.gov/RefSeq/key.html

#accessions
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
	     [2157, 9], #archaea, get WGS
	     [2,    9], #bacteria, get WGS
	     [2759, 9], #euks, not ready for WGS -- lots of data with minimal annotations.  Will get there though!
	     [10239,9], #viruses, phages,
	     [12884,9], #viroids
	    );

my $pm = new Parallel::ForkManager($forks);
foreach my $item (@taxids)
  {
    my ($id, $type) = @$item;
    my $url = $base_url."type=$type"."&taxid=$id";
    print "Fetching NCBI Genomes List from HTTP: $url\n";
    my $content =get($url);
    my @tables = split /<\/table>/, $content;
    my @rows = split /<tr>/,$tables[8];
    my %data;
    foreach my $row (@rows)
      {

 	next unless $row =~/^<td>/;
 	$row =~ s/\n|\r//g;
 	my @cols = split /<td.*?>/, $row;
 	foreach my $item (@cols)
 	  {
 	    $item =~ s/<.*?>//g;
 	  }
 	next unless $cols[3];
 	next if $cols[3] =~ /accession/;
 	$cols[2] =~ s/&nbsp;//g;
 	$cols[2] =~ s/chromosome//i;
 	$cols[2] =~ s/chr//i;
 	$cols[2] =~ s/^\s+//;
 	$cols[2] =~ s/\s+$//;
 	$cols[2] = undef if $cols[2] && ($cols[2] eq " " || $cols[2] =~ /nbsp/);
 	$cols[2] = "unknown" unless $cols[2];
 	#$data{organism}{chromosome} = [accns]
 	push @{$data{$cols[1]}{$cols[2]}}, $cols[3];

      }
#    print Dumper \%data;
    my $orgs = process_orgs(\%data);
#    print Dumper $orgs;
    foreach my $org (keys %$orgs)
      {
	foreach my $proj (keys %{$orgs->{$org}})
	  {
	    my @accns = map {$_->{accn}} @{$orgs->{$org}{$proj}};
	    process_genome_and_load(genome=>{orgs=>[$org],accns=>[@accns]});
	  }
      }
  }
#my $genomes = get_NCBI_genomes();

if (@skipped)
  {
    print "Commands that had automatially skipped datasets:\n";
    print join ("\n", @skipped),"\n";
  }

exit();

sub process_genome_and_load
  {
    my %opts = @_;
    my $genome = $opts{genome};
    my $prog = '/home/elyons/projects/CoGeX/scripts/load_genomes_n_stuff/genbank_genome_loader.pl';
    $prog .= " -autoupdate" if $autoupdate;
    #need to break apart organelle genomes from nuclear genomes
    my @tmp; #temp place to store genomes in case some need to be processed separately in addition to the whole genome.  Such as is the case with mitochondria, chloroplasts, plastids.
    if (@{$genome->{orgs}} > 1) #don't need to check if there is only one genome for the organism
      {
	my @subgenome_names = qw(mitochondrion mitochondrial chloroplast apicoplast plastid);
	my $i=0; #array index;
	my %no_organelles= (org=>[], accns=>[]); #storage for genomes without organelles
	foreach my $org (@{$genome->{orgs}})
	  {
	    my $organelle_found=0;
	    foreach my $check (@subgenome_names, "whole genome shotgun sequencing project")
	      {
		if ($org =~ /$check/)
		  {
		    $organelle_found=1;
		    #add organelle as separate genome
		    push @tmp, {
				org=>[$org],
				accns=>[$genome->{accns}[$i]],
			       };
		  }
	      }
	    unless ($organelle_found)
	      {
		push @{$no_organelles{org}}, $org;
		push @{$no_organelles{accns}}, $genome->{accns}[$i];
	      }
	    $i++;
	  }
	push @tmp, \%no_organelles;
      }
    else #only one genome
      {
	push @tmp, $genome;
      }

    foreach my $genome (@tmp)
      {
	#	$pm->start and next;
	my $run = $prog." -accn ".join (" -accn ", @{$genome->{accns}});
	#$run .= " -chr '$cols[2]'" if $cols[2];
	$run .= " -td '/tmp/gb/'";
	$run .= " -go";
	$run .= " -autoskip";
	$run .= " -delete_src_file";
	print "\n",$run,"\n";
	my $previously_loaded = 0;
	my $skipped = 0;
#	next;
	open (IN, $run." |");
	while (<IN>)
	  {
	    $previously_loaded =1 if /previously loaded/;
	    $skipped=1 if /skipping/;
	    print $_;
	  }
	close IN;
	push @skipped, $run if $skipped;
	if ($pause && !$previously_loaded)
	  {
	    print "Paused:  press enter to continue\n";
	    <STDIN>;
	  }
	#	$pm->finish;
      }
  }
sub process_orgs
  {
    my $data = shift;
    my %orgs;
    foreach my $org(keys %$data)
      {
	my $check = 0;
	foreach my $chr (keys %{$data->{$org}})
	  {
	    if (scalar @{$data->{$org}{$chr}} > 1 || $check)
	      {
#		print Dumper $data->{$org}{$chr},"\n";
		foreach my $accn (@{$data->{$org}{$chr}})
		  {
		    my $project_id_hash = get_project($accn);
		    my $project_id = $project_id_hash->{$accn};
		    next unless $project_id;
		    push @{$orgs{$org}{$project_id}}, {accn=>$accn};
		  }
		$check =1;
	      }
	    else
	      {
		push @{$orgs{$org}{"NA"}}, {accn=>$data->{$org}{$chr}[0]};
	      }
	  }
      }
    return \%orgs;
  }

sub get_project  {
    my $accn = shift;
    $accn = join ",", @$accn if ref ($accn) =~ /array/i;
    my $summary_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=text&complexity=3&id=";
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
#	next;
      }

    my %data;
#    while ($summary =~/ACCESSION\s+(\S+).*?DBLINK\s+Project:(\d+)/gxs)
#    print $summary if $accn eq "AC_000167";
    while ($summary =~ /DBLINK\s+Project:(\d+)/g)
      {
	$data{$accn} = $1;
      }
#    print scalar keys %data,"\n";
    return \%data;
  }

sub get_NCBI_genomes
  {
    my $esearchgenomeprj = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genomeprj&term=all%5Bfilter%5D&retmax=999999"; #genomeprj let's me assemble genomes when there are multiple assemblies (no viruses)
    my $esearchgenome = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=all%5Bfilter%5D&retmax=999999"; #get all genome ids
    my @genomes;

    #get genome projects
    my $entry = get($esearchgenomeprj);
    my %projects; #store genome_id=genomeprj_id
    my $i =0;
    while ($entry =~ /<Id>(\d+)<\/Id>/ig)
      {
#	last if $i == 2000;
	my $gpid = $1;
	print "working on $gpid\n";
	print get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=genomeprj&retmode=text&complexity=0&id=$gpid");

	my $gids = get_genomes_for_genomeprj($gpid); #metagenomes, etc won't have genomes. . .
	my @gids;
	foreach my $gid (@$gids)
	  {
	    next unless $gid; #many genome projects have no specific genomes -- metagenomes, ENCODE
	    $projects{$gid}=$gpid;
	    push @gids, $gid;
	  }
	if (@gids)
	  {
#	    next; #this is temp to test other routine
	    my $gids = join ",", @gids;
	    my ($accns, $orgs) = get_genome_accns($gids);
	    my $res = {gids=>\@gids,
		       gpid=>$gpid,
		       accns=>$accns,
		       orgs=>$orgs};
	    process_genome_and_load(genome=>$res);
#	    push @genomes, $res;
	  }
	else
	  {
#	    next;
	    $gpid = 38683;
	    my $ncids = get_nuccore_for_genome_project($gpid);
	    unless (@$ncids)
	      {
		print "\tunable to get genome or nuccore for genome project id $gpid\n";
		next;
		#trying to figure someway out to get the data.
		get_all_genome_prj_links($gpid);
		foreach my $id (@$ncids)
		  {
		    print "\tfetching nuccore for genome project $id\n";
		    get_nuccore_for_genome_project($id);
		  }
	      }
#	    my %accns;
	    my %orgs;
	    foreach my $ncid (@$ncids)
	      {
		my ($accns) = get_accns_for_nuccore($ncid);
		my ($org) = get_organism_for_nuccore($ncid);
		foreach my $accn (@$accns)
		  {
		    $orgs{$org}{$accn}=1;
		  }
	      }
	    foreach my $org (keys %orgs)
	      {
		my $res = {gids=>\@gids,
			   gpid=>$gpid,
			   ncids=>$ncids,
			   accns=>[keys %{$orgs{$org}}],
			   orgs=>[$org]
			  };

		print Dumper $res;
		<STDIN>;
		process_genome_and_load(genome=>$res);
	      }
#	    push @genomes, $res;
#	    print "\n";
	  }
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

#	my $gids = join ",", @{$genome->{gids}};
#	my ($accns, $orgs) = get_genome_accns($gids);
#	$genome->{accns}=$accns;
#	$genome->{orgs}=$orgs;# if @$accns;
	$i++;
      }
    return \@genomes;
  }

sub get_all_genome_prj_links
    {
      my $id = shift;
      my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=all&dbfrom=genomeprj&retmode=text&id=";
      my $entry = get ($url.$id);
#      print $entry;
    }

sub get_genomes_for_genomeprj
    {
      my $gid = shift;
      my $elink = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=genome&dbfrom=genomeprj&id=$gid";
      print "running get_genomes_for_genomeprj: $elink\n";
      my $entry = get($elink);
      my @ids;
      return \@ids unless $entry;
      my $xml = XMLin($entry);
#      print $entry;
      my $items = ref ($xml->{LinkSet}{LinkSetDb}{Link}) =~ /array/i ? [@{$xml->{LinkSet}{LinkSetDb}{Link}}] : [$xml->{LinkSet}{LinkSetDb}{Link}];

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
      print "running get_genomes_accns\n";
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
	  push @orgs, $1;
	}
#      print "successfully got accessions:\n";
#      print "\t", join ("\t", @orgs),"\n";
#      print "\t", join ("\t", @accns),"\n";
      return \@accns, \@orgs;
    }

sub get_nuccore_for_genome_project
   {
     my $gpid = shift;
     my $elink = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=nuccore&dbfrom=genomeprj&id=$gpid";
     print "running get_nuccore_for_genome_project: $elink\n";
     my $entry = get($elink);
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

sub get_accns_for_nuccore
   {
     my $id = shift;

     my $esumg = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=text&complexity=0&id=$id";
     print "running get_accns_for_nuccor: $esumg\n";
     my $item = get($esumg);
#     print $item;
     my @accns;
     my @orgs;
     while ($item=~ /<Item Name="Caption" Type="String">(.*?)<\/Item>/g)
       {
	 push @accns, $1;
       }
     while ($item=~ /<Item Name="Title" Type="String">(.*?)<\/Item>/g)
       {
#	 print $1,"\n";
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

sub get_organism_for_nuccore
   {
     my $nc = shift;
     my $elink = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&id=$nc";
     print "running get_nuccore_for_genome_project: $elink\n";
     my $entry = get($elink);
     return unless $entry;
     my $xml = XMLin($entry);
     my $items = ref ($xml->{LinkSet}{LinkSetDb}{Link}) =~ /array/i ? [@{$xml->{LinkSet}{LinkSetDb}{Link}}] : [$xml->{LinkSet}{LinkSetDb}{Link}];
     my $esummary = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&retmode=text&complexity=0&id=".$items->[0]{Id};
     print "running get_nuccore_for_genome_project: $esummary\n";
     my $entry2 = get($esummary);
     my $org_name;
     ($org_name) = $entry2=~ /<Item Name="ScientificName" Type="String">(.*?)<\/Item>/;
     return $org_name;
   }
