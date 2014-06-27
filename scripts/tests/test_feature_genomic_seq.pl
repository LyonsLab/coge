#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my ($help, @orgids, $version, $type, $chr, $DEBUG, $length, $org_search);

GetOptions ("h|help" =>  \$help,
            "o|org=s" => \@orgids,
            "v|version=s" => \$version,
            "t|type=s"    => \$type,
            "d|debug"     => \$DEBUG,
            "c|chr=s"     => \$chr,
	    "chr_length|cl=s" => \$length,
	    "org_search|os=s" =>\$org_search,
            );

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
$type = "CDS" unless $type;

my @aa = qw(R H K D E S T N Q C G P A I L M F W Y V);
my @codons = gen_codon_table();
print STDERR "length $length\n" if $DEBUG;

$| =1;
print join ("\t", qw (GROUP TYPE NAME LENGTH CDS_COUNT SYS_1 SYS_2), (map {$_, $_."%"} sort @aa), (map {$_, $_."%"} sort @codons), (map {$_, $_."%"} qw(GC AT) ) ),"\n";
my @orgs;
unless (@orgids)
  {
    @orgs = $coge->resultset('Organism')->all;
  }
else
  {
    foreach my $orgid (@orgids)
      {
	foreach my $org ($coge->resultset('Organism')->find($orgid))
	  {
	    push @orgs, $org;
	  }
      }
  }

foreach my $org (@orgs)
  {
    if ($org_search) {next unless $org->name =~ /$org_search/i || $org->description =~ /$org_search/i};
    my $oname = $org->name." ".$org->description;
#    print $oname,"\n";
#    next;

    my $otype = get_type($oname);
    my $group = get_group($oname);
#    my $search ={organism_id=>$org->id};
#    $search->{version}=$version if $version;
    foreach my $ds ($coge->get_current_datasets_for_org(org=>$org->id))
      {
	print STDERR "working on ",$ds->name,"\n" if $DEBUG;
	foreach my $feat ($ds->features)
	  {
        my $seq1 = $feat->genomic_sequence();
        my $seq2 = $feat->genomic_sequence_old();
        print "UNEQUAL\n" unless $seq1 eq $seq2;
        print "PASSED\n" if $seq1 eq $seq2;
	  }
      }
 }

sub get_group
  {
    my $name = shift;
    my $group;
    if ($name =~ /virus/i)
      {
	$group = "Virus";
      }
    elsif ($name =~ /mitochondr/i)
      {
	$group = "mitochondrion";
      }
    elsif ($name =~ /chloroplast/i)
      {
	$group = "chloroplast";
      }
    elsif ($name =~ /Archaea/i)
      {
	$group = "Archaea";
      }
    elsif ($name =~ /Bacteria/i)
      {
	$group = "Bacteria";
      }
    elsif ($name =~ /Eukary/i)
      {
	$group = "Eukaryote";
      }
    else
      {
	$group = "unknown";
      }
    return $group;
  }

sub get_type
  {
    my $name = shift;
    my $virus_type;
    if ($name =~ /dsDNA/i)
      {
	$virus_type = "dsDNA";
      }
    elsif ($name =~ /ssDNA/i)
      {
	$virus_type = "ssDNA";
      }
    elsif ($name =~ /ssRNA/i)
      {
	$virus_type = "ssRNA";
      }
    elsif ($name =~ /dsRNA/i)
      {
	$virus_type = "dsRNA";
      }
    elsif ($name =~ /retrovir/i)
      {
	$virus_type = "RT RNA: retrovirus";
      }
    elsif ($name =~ /hepadna/i)
      {
	$virus_type = "RT DNA: hepadna";
      }
    elsif ($name =~ /single stranded DNA/i)
      {
	$virus_type = "ssDNA";
      }
    elsif ($name =~ /single stranded RNA/i)
      {
	$virus_type = "ssRNA";
      }
    elsif ($name =~ /Caulimo/i)
      {
	$virus_type = "RT DNA: caulimo";
      }
    elsif ($name =~ /RNA/i)
      {
	$virus_type = "RNA";
      }
    elsif ($name =~ /DNA/i)
      {
	$virus_type = "DNA";
      }
    else
      {
	$virus_type = "unknown";
      }
    return $virus_type;
  }

sub gen_codon_table
  {
    my @nt = qw (A T C G);
    my @codons;
    foreach my $nt1 (sort @nt)
      {
	foreach my $nt2 (sort @nt)
	  {
	    foreach my $nt3 (sort @nt)
	      {
		push @codons, $nt1.$nt2.$nt3;
	      }
	  }
      }
    return @codons;
  }
