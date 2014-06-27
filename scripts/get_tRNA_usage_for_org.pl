#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use Data::Dumper;

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
$type = "tRNA" unless $type;
my $aragorn = "";

my @aa = qw(R H K D E S T N Q C G P A I L M F W Y V);
my %codons = map {$_,0} gen_codon_table();

print STDERR "length $length\n" if $DEBUG;

$| =1;
print join ("\t", qw (GROUP TYPE NAME), (map {$_, $_."%"} sort keys %codons)),"\n";
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
    my $org_length = 0;
    my %codon_usage = %codons;;
    my $codon_total = 0;
    foreach my $ds ($coge->get_current_datasets_for_org(org=>$org->id))
      {
	print STDERR "working on ",$ds->name,"\n" if $DEBUG;
#	print "working on ",$ds->name,"\n";
	if ($length)
	  {
	    my $total_length =0;
	    foreach my $chr ($ds->get_chromosomes)
	      {
		my $last = $ds->last_chromosome_position($chr);
		$total_length += $last;
		$org_length += $last;
	      }
	    print STDERR "Checking length:  limit: $length, total length $total_length\n"  if $DEBUG;
	    next unless $total_length;
	    next if $total_length > $length;
	  }
	my $file = gen_fasta($ds);
	my $data = run_aragorn($file);
	foreach my $c (keys %codon_usage)
	  {
	    next unless $data->{$c};
	    $codon_usage{$c}+=$data->{$c};
	    $codon_total+=$data->{$c};
	  }
#	foreach my $feat ($coge->resultset('Feature')->search({dataset_id=>$ds->id, "feature_type.name"=>{like=>"%trna%"}},{join =>"feature_type"}))
# 	  {
# 	    my $anti_codon;
# 	    my $codon;
# 	    my $aa;
# 	    my $annotations = join ("\t",map {$_->annotation} $feat->annotations);
# 	    ($anti_codon) = $annotations =~ /anticodon:?\s*(\w\w\w)/;
# 	    ($codon) = $annotations =~ /codon recognized:?\s*(\w\w\w)/;
# 	    if ($codon)
# 	      {
# 		$anti_codon = $codon;
# 		$anti_codon =~ tr/atcgATCG/tagcTAGC/;
# 	      }
# 	    $codon_usage{$anti_codon}++ if $anti_codon;
# 	    $anti_codon = $annotations unless $anti_codon;
# 	    print ">".join ("\t", $group, $otype, $oname, $feat->type->name,$anti_codon),"\n";
# #	    print $feat->genomic_sequence,"\n";
# 	  }
#       }
#     foreach (@codons)
#       {
# 	$codon_usage{$_} = 0 unless $codon_usage{$_};
#       }
#     print join ("\n", map {$_ ." " .$codon_usage{$_}} @codons),"\n";
#     next;
#     my $codon_total =0;
#     map {$codon_total+=$codon_usage{$_}} sort @codons;
      }
    print join ("\t", $group, $otype, $oname, (map {$codon_usage{$_}, sprintf("%.4f",$codon_usage{$_}/$codon_total)} sort keys %codons)),"\n" if $codon_total;
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

sub gen_fasta
  {
    my $ds = shift;
#    my $title = join (", ", join (", ", map {"chr:".$_." v:".$ds->version." ds:".$ds->id} $ds->get_chromosomes));
#    my $md5 = md5_hex($title);
    my $file = "tmp.fasta";
    open (OUT, ">$file") || die "Can't open $file for writing: $!";;
    foreach my $chr (sort $ds->get_chromosomes)
      {
	my $title =  $ds->organism->name." (v". $ds->version.") "."chromosome: $chr".", CoGe database id: ".$ds->id;
	$title =~ s/^>+/>/;
	print OUT ">".$title."\n";
	print OUT $ds->get_genomic_sequence(chr=>$chr),"\n";
      }
    close OUT;
    return $file if -r $file;
    return 0;
  }

sub run_aragorn
  {
    my $file = shift;
    my $cmd = $aragorn ." -t -w ".$file;
    my %data;
    print STDERR "running $cmd" if $DEBUG;
    open (IN, $cmd."|");

    while (<IN>)
      {
	next if /^>/;
#	chomp;
#	my @line = split/\s+/;
	my ($anticodon) = /\((.*)\)/;
	next unless $anticodon;
#	print $_,$anticodon,"\n";
       $anticodon =~ tr /atcgATCG/tagcTAGC/;
	$data{uc($anticodon)}++;

      }
    close IN;
    return \%data;
  }
