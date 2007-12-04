#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use DBIxProfiler;

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

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

$type = "CDS" unless $type;
my $aragorn = "/home/elyons/bin//aragorn";



my @aa = qw(R H K D E S T N Q C G P A I L M F W Y V);
my @codons = gen_codon_table();
print STDERR "length limit $length\n" if $DEBUG && $length;

$| =1;
print join ("\t", qw (GROUP TYPE NAME LENGTH CDS_COUNT SYS_1 SYS_2), (map {$_, $_."%"} sort @aa), (map {$_, $_."%"} sort @codons), (map {"t".$_, "t".$_."%"} sort @codons), (map {$_, $_."%"} qw(GC AT) ) ),"\n";
my @orgs;
unless (@orgids)
  {
    @orgs = $coge->resultset('Organism')->all;
  }
else
  {
    foreach my $orgid (@orgids)
      {
	foreach my $org ($coge->resultset('Organism')->resolve($orgid))
	  {
	    print $org->name, " matched $orgid\n" if $DEBUG;
	    push @orgs, $org;
	  }
      }
  }

foreach my $org (@orgs)
  {
    if ($org_search) {next unless $org->name =~ /$org_search/i || $org->description =~ /$org_search/i};
    my $oname = $org->name;
    $oname .= ": ".$org->description if $org->description;
#    print $oname,"\n";
#    next;

    my $otype = get_type($oname);
    my $group = get_group($oname);
#    my $search ={organism_id=>$org->id};
#    $search->{version}=$version if $version;
    my ($sys1, $sys2) = (0,0);
    my $cds_count = 0;
    my $org_length = 0;
    my %aa_usage = map {$_,0} @aa;
    my %codon_usage = map {$_, 0} @codons;
    my %trna_codon_usage = map {$_,0} @codons;
    my $trna_count = 0;
    my $at = 0;
    my $gc = 0;
#    foreach my $ds ($coge->resultset('Dataset')->search($search))
    foreach my $ds ($coge->get_current_datasets_for_org(org=>$org->id))
      {
#	next unless $ds->id == "6002";
	print STDERR "working on ",$oname,": ",$ds->name,"\n" if $DEBUG;
	my $total_length =0;
	foreach my $chr ($ds->get_chromosomes)
	  {
	    my $last = $ds->last_chromosome_position($chr);
	    $total_length += $last;
	    $org_length += $last;
	  }
	print STDERR "Checking length:  limit: $length, total length $total_length\n"  if $DEBUG && $length;
	next unless $total_length;
	next if $length && $total_length > $length;
	#tRNA stuff
	my $file = gen_fasta($ds);
	my $trna_data = run_aragorn($file);
#	print join ("\n", map{$_."=>".$trna_data->{$_}} sort keys %$trna_data);
	foreach my $c (keys %trna_codon_usage)
	  {
	    next unless $trna_data->{$c};
	    $trna_codon_usage{$c}+=$trna_data->{$c};
	    $trna_count+=$trna_data->{$c};
	  }
	#feature specific stuff
	print STDERR "processing features: " if $DEBUG;
	foreach my $feat ($ds->features({"feature_type.name"=>$type},{join=>"feature_type"}))
	  {
	    next unless $feat->id =~ /\d+/;
#	    print join ("\t", $feat->id, $feat->names, $feat->type->name),"\n";
#	    print Dumper $feat->id;
#	    next if $type && $feat->type->name ne $type;
	    $cds_count++;
	    #synthase system
	    my ($code1, $code2) = $feat->percent_translation_system(counts=>1);
	    $sys1 += $code1;
	    $sys2 += $code2;
	    #aa usage
	    my $aa_usage = $feat->aa_frequency(counts=>1);
	    foreach my $aa (keys %$aa_usage)
	      {
		$aa_usage{uc($aa)}+=$aa_usage->{$aa};
	      }
	    #codon usage
	    my ($codon_usage) = $feat->codon_frequency(counts=>1);
	    foreach my $codon (keys %$codon_usage)
	      {
		$codon_usage{uc($codon)}+=$codon_usage->{$codon};
	      }
	    my @gc_at = $feat->gc_content(counts=>1);
	    unless ($gc_at[0] || $gc_at[1])
	      {
		print $feat->id, " had no sequence\n";
	      }
	    $gc+= $gc_at[0] if $gc_at[0] && $gc_at[0] =~ /^\d+$/;
	    $at+= $gc_at[1] if $gc_at[1] && $gc_at[1] =~ /^\d+$/;
	    print STDERR "." if ($cds_count%10 == 0 && $DEBUG);
	  }
	print STDERR "finished!\n" if $DEBUG;
      }
    
    if (($gc+$at) == 0)
      {
	print STDERR $oname ." has no sequence\n" if $DEBUG;
	next;
      }

    my $aa_total =0;
    map {$aa_total+=$aa_usage{$_}} sort @aa;
    my $codon_total =0;
    map {$codon_total+=$codon_usage{$_}} sort @codons;
    my $nt_total = $at+$gc;
    print join ("\t", $group, $otype, $oname, $org_length, $cds_count);
    print "\t";
    if ($sys1 || $sys2)
      {
	print join ("\t", map {sprintf("%.4f", $_)} $sys1/($sys1+$sys2), $sys2/($sys1+$sys2));
      }
    else
      {
	print "\tN/A\tN/A";
      }
    if ($aa_total)
      {
	print "\t",join ("\t", (map {$aa_usage{$_}, sprintf("%.4f",$aa_usage{$_}/$aa_total)} sort @aa));
      }
    else
      {
	print "\tN/A\tN/A"x scalar@aa;
      }
    if ($codon_total)
      {
	print join "\t",("\t", (map {$codon_usage{$_}, sprintf("%.4f",$codon_usage{$_}/$codon_total)} sort @codons));
      }
    else
      {
	print "\tN/A\tN/A"x scalar@codons;
      }
    if ($trna_count)
      {
	print "\t",join ("\t", (map {$trna_codon_usage{$_}, sprintf("%.4f",$trna_codon_usage{$_}/$trna_count)} sort @codons));
      }
    else
      {
	print "\t0\t0"x scalar@codons;
      }
    print "\t",join ("\t", (map {$_, sprintf("%.4f",$_/$nt_total)} $gc, $at)),"\n";
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
	$group = "Mitochondrion";
      }
    elsif ($name =~ /chloroplast/i)
      {
	$group = "Chloroplast";
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
    print STDERR "running $cmd\n" if $DEBUG;
    my $count=0;
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
	$data{reverse(uc($anticodon))}++;
	$count++;
      }
    print STDERR "\tfound $count tRNAs\n";
    close IN;
    return \%data;
  }
    
