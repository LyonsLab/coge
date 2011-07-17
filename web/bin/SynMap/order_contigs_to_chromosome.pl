#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGeX;
use Text::Wrap;
use CGI;
use CoGe::Accessory::Web;
use CoGe::Accessory::SynMap_report;

use vars qw($synfile $coge $DEBUG $join $FORM $P $GZIP $GUNZIP);


$P = CoGe::Accessory::Web::get_defaults();
$GZIP = $P->{GZIP};
$GUNZIP = $P->{GUNZIP};
GetOptions (
	    "debug"=>\$DEBUG,
	    "file|f=s"=>\$synfile, # file_name.aligncoords from SynMap
	    "join=i"=>\$join, #is the output sequence going to be joined together using "N"s for gaps.  Set to a number to be true, and whatever number it is will be the number of N's used to join sequence.  
	   );
$FORM = new CGI;
$synfile = "/opt/apache/".$FORM->param('f') if $FORM->param('f');

$join = 100 unless defined $join;
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

my $synmap_report = new CoGe::Accessory::SynMap_report;
$synfile = gunzip($synfile);
my ($chr1, $chr2) = $synmap_report->parse_syn_blocks(file=>$synfile);
gzip($synfile);
#print $FORM->header;
print qq{Content-Type: text/plain

};
#print "<pre>";
process_sequence($chr1);
#print "</pre>";

sub process_sequence
  {
    my $chrs = shift;
    my %dsg; #store CoGe dataset group objects so we don't have to create them multiple times
    my $count = 0; #number of blocks processed per matched chromosome
    my %out_chrs; #seen chromosomes for printing, let's me know when to start a new fasta sequence
    my %in_chrs; #seen chromosomes coming in, need to use this to identify those pieces that weren't used and to be lumped under "unknown"
    my $seq; #sequence to process and dump
    my $header; #header for sequence;
    foreach my $item (@$chrs)
      {
	my $chr = $item->{chr};
	my $out_chr = $item->{match};
	my $dsgid = $item->{dsgid};
	my $dsg = $dsg{$dsgid};
	$dsg = $coge->resultset('DatasetGroup')->find($dsgid) unless $dsg;
	$dsg{$dsg->id} = $dsg;
	my $strand = $item->{rev} ? -1 : 1;
	if ($seq && !$out_chrs{$out_chr})
	  {
	    #we have a new chromosome.  Dump the old and get ready for the new;
	    print_sequence(header=>$header, seq=>$seq);
	    $count=0;
	    $seq=undef;
	  }
	if ($join)
	  {
	    if ($count)
	      {
		$seq .= "N"x$join;
	      }
	    else#need to print fasta header
	      {
		$header = ">$out_chr";#use the matched organism's chromosome
#		$header .= $dsg->organism->name;
#		$header .= " ".$dsg->organism->description if $dsg->organism->description;
#		$header .= " ".$dsg->name if $dsg->name;
#		$header .= ": ".$dsg->description if $dsg->description;
#		$header .= " (v".$dsg->version." ".$dsg->type->name.")\n";
	      }
	    $seq .= $dsg->genomic_sequence(chr=>$chr, strand=>$strand);
	    $count++;
	    $out_chrs{$out_chr}++;
	    $in_chrs{uc($chr)}++;
	  }
	else
	  {
	    print $dsg->fasta(chr=>$chr);
	  }
      }
    if ($seq)
      {
	print_sequence(header=>$header,seq=>$seq);
      }
    #need to get all the pieces that didn't fit
    $header = ">Unknown";
    $seq = undef;
    $count=0;
    foreach my $dsg (values %dsg)
      {
	foreach my $chr ($dsg->chromosomes)
	  {
	    next if $in_chrs{uc($chr)};
	    $seq .= "N"x$join if $count;
	    $seq .= $dsg->genomic_sequence(chr=>$chr);
	    $count++;
	  }
      }
    print_sequence(header=>$header, seq=>$seq);
  }


sub print_sequence
    {
      my %opts = @_;
      my $header = $opts{header};
      my $seq = $opts{seq};
      $Text::Wrap::columns=80;
      print $header,"\n";
      print wrap('','',$seq),"\n";
    }

sub gzip
    {
      my $file = shift;
      return $file unless $file;
      return $file if $file =~ /\.gz$/;
      return $file.".gz" if -r "$file.gz";
      `$GZIP $file` if -r $file;
      my $tmp = $file.".gz";
      return -r $tmp ? $tmp : $file;
    }

sub gunzip
    {
      my $file = shift;
      return $file unless $file;
      return $file unless $file =~ /\.gz$/;
      `$GUNZIP $file` if -r $file;
      my $tmp = $file;
      $tmp =~ s/\.gz$//;
      return -r $tmp ? $tmp : $file;
    }
