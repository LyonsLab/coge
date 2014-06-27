#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use DBIxProfiler;
use DBI;

my ($help, @orgids, $version, $type, $chr, $DEBUG, $length, $org_search, $sqlite, $fasta, $short);

GetOptions ("h|help" =>  \$help,
            "o|org=s" => \@orgids,
            "v|version=s" => \$version,
            "t|type=s"    => \$type,
            "d|debug"     => \$DEBUG,
            "c|chr=s"     => \$chr,
	    "org_search|os=s" =>\$org_search,
	    "fasta"       => \$fasta,
	    "short"=>\$short,
            );

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);
$type = "CDS" unless $type;
my $aragorn = "/home/elyons/bin//aragorn";

$| =1;

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

unless (@orgs)
  {
    print "No organisms identified\n";
    exit;
  }
print join ("\t", qw{ORG DSID CHR START STOP STRAND AA ANTICODON CODON SEQ}),"\n" unless $fasta;

foreach my $org (@orgs)
  {
    if ($org_search) {next unless $org->name =~ /$org_search/i || $org->description =~ /$org_search/i};
    my $oname = $org->name;
    $oname .= ": ".$org->description if $org->description;
    foreach my $ds ($org->current_datasets())
      {
	#tRNA stuff
	my $file = gen_fasta($ds);
	my $trna_data = run_aragorn($file);
	if ($fasta)
	  {
	    foreach my $org(sort keys %$trna_data)
	      {
		foreach my $aa (sort keys %{$trna_data->{$org}})
		  {
		    my $count =1;
		    foreach my $seq (sort keys %{$trna_data->{$org}{$aa}})
		      {
			my $head;
			foreach my $item (split/\s+/,$org,2)
			  {
			    $head.=substr($item,0,1);
			  }
			$head .= $trna_data->{$org}{$aa}{$seq}{aa};
			$head .= $count;
			$head .= $trna_data->{$org}{$aa}{$seq}{codon};
			print ">$head\n";
			print $seq,"\n";
			$count++;
		      }
		  }
	      }
	  }
      }
  }
sub gen_fasta
  {
    my $ds = shift;
    my $file = "tmp.fasta";
    open (OUT, ">$file") || die "Can't open $file for writing: $!";;
    foreach my $chr (sort $ds->get_chromosomes)
      {
	open (OUT, ">$file");
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
    my $org;
    my $chr;
    my $dsid;
    my $ds;
    while (<IN>)
      {
#	print $_;
	chomp;
	if (/^>/)
	  {
	    $org = $_;
	    $org =~ s/^>//;
	    ($dsid) = /id:\s(\d+)/;
	    ($chr) = /chromosome: (\w+),/;
	    $ds = $coge->resultset('Dataset')->find($dsid);
	    next;
	  }
	my @line = split/\s+/;
	my ($aa,$anticodon) = /(\w+)\((\w+)\)/;
	my ($start, $stop) = $line[1]=~/(\d+),(\d+)/;
	my $strand = $line[1]=~/c/ ? -1 : 1;
	next unless $anticodon;
#	print $_,$anticodon,"\n";
        $anticodon =~ tr /atcgATCG/tagcTAGC/;
	my $codon = reverse(uc($anticodon));
	$count++;
	my $seq = $ds->genomic_sequence(start=>$start,
					stop=>$stop,
					strand=>$strand,
					chr=>$chr,
					);
	if ($data{$org}{$aa}{$seq})
	  {
	    $data{$org}{$aa}{$seq}{count}++;
	  }
	else
	  {
	    $data{$org}{$aa}{$seq}={
			 count=>1,
			 org=>$org,
			 dsid=>$dsid,
			 chr=>$chr,
			 start=>$start,
			 stop=>$stop,
			 strand=>$strand,
			 aa=>$aa,
			 anticodon=>$anticodon,
			 codon=>$codon,
			};
	  }
	unless ($fasta)
	  {
	    print join ("\t",$org, $dsid, $chr, $start, $stop, $strand, $aa, $anticodon, $codon, $seq),"\n";
	  }

      }
    print STDERR "\tfound $count tRNAs\n";
    close IN;
    return \%data;
  }
