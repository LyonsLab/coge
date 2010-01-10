#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGeX;
use Text::Wrap;
use CGI;

use vars qw($synfile $coge $DEBUG $join $FORM);

GetOptions (
	    "debug"=>\$DEBUG,
	    "file|f=s"=>\$synfile, # file_name.aligncoords from SynMap
	    "join=i"=>\$join, #is the output sequence going to be joined together using "N"s for gaps.  Set to a number to be true, and whatever number it is will be the number of N's used to join sequence.  
	   );
$FORM = new CGI;
$synfile = "/opt/apache/".$FORM->param('f') if $FORM->param('f');

$join = 100 unless defined $join;
$coge = CoGeX->dbconnect();

my ($ordered_chrs) = parse_syn_blocks($synfile);
print $FORM->header;
print "<pre>";
print_sequence($ordered_chrs);
print "</pre>";

sub print_sequence
  {
    my $chrs = shift;
    my %dsg;
    my $count = 0;
    my $seq;
    foreach my $item (@$chrs)
      {
	my ($dsgid, $chr) = split/_/, $item->{chr};
	$dsgid =~ s/^\D//;
	my $dsg = $dsg{$dsgid};
	$dsg = $coge->resultset('DatasetGroup')->find($dsgid) unless $dsg;
	$dsg{$dsg->id} = $dsg;
	my $strand = $item->{rev} ? -1 : 1;
	if ($join)
	  {
	    if ($count)
	      {
		$seq .= "N"x$join;
	      }
	    else#need to print fasta header
	      {
		print ">";
		print $dsg->organism->name;
		print " ".$dsg->organism->description if $dsg->organism->description;
		print " ".$dsg->name if $dsg->name;
		print ": ".$dsg->description if $dsg->description;
		print " (v".$dsg->version." ".$dsg->type->name.")\n";
	      }
	    $seq .= $dsg->genomic_sequence(chr=>$chr, strand=>$strand);
	  }
	else
	  {
	    print $dsg->fasta(chr=>$chr);
	  }
	$count++;
      }
    if ($seq)
      {
	$Text::Wrap::columns=80;
	print wrap('','',$seq);
      }
  }

sub parse_syn_blocks
  {
    my $file = shift;


    my $blocks1=[];
    my $blocks2=[];
    open (IN, $file) || die "Can't open $file for reading: $!";
    $/ = "\n#";
    while (<IN>) #get blocks
      {
	next unless $_;
	s/#//g;
	my ($block1, $block2) = process_syn_block($_);
	push @$blocks1 , $block1;
	push @$blocks2 , $block2;
      }
    close IN;
    $/="\n";
    #which organism has more pieces (aka chromosomes, contigs, etc.)  we are going to assemble the one with more pieces
#    print Dumper \@blocks;
    my $chrs1={};
    my $chrs2={};
    foreach my $item (@$blocks1)
      {
	$chrs1->{$item->{name}}{count}++;
	$chrs1->{$item->{name}}{score}+=$item->{score};
      }
    foreach my $item (@$blocks2)
      {
	$chrs2->{$item->{name}}{count}++;
	$chrs2->{$item->{name}}{score}+=$item->{score};
      }
#    map {$chrs1->{$_->{name}}{count}++} @$blocks1;
#    map {$chrs2->{$_->{name}}{count}++} @$blocks2;
    #blocks1 will contain fewer chromosomes; blocks2 will be ordered by it.
    my $switched =keys %$chrs1 > keys %$chrs2 ? 1 : 0;
    if ($switched)
      {
	($blocks1, $blocks2) = ($blocks2, $blocks1);
	($chrs1, $chrs2) = ($chrs2, $chrs1);
      }

    my $ordered1 =[]; #storage for ordered chromosomes
    my $ordered2 =[]; #storage for ordered chromosomes
    my %seen;
    foreach my $chr1 (sort{$chrs1->{$b}{score} <=> $chrs1->{$a}{score}} keys %$chrs1)
      {
	push @$ordered1, {chr=>$chr1};
	my @blocks;
	for (my $i=0; $i< @$blocks1; $i++)
	  {
	    my ($block1, $block2) = ($blocks1->[$i],$blocks2->[$i]);
	    next if $block1->{name} ne $chr1;
	    push @blocks, $block2;
	  }
	#Need to check if a given chromosome in @blocks occurs more than once.  This will happen if there is a segmental duplication, or something along those lines.
	my %block_check;
	foreach my $block (@blocks)
	  {
	    if ($block_check{$block->{name}})
	      {
		$block_check{$block->{name}}=$block if $block->{score} > $block_check{$block->{name}}{score}; #higher score wins
	      }
	    else
	      {
		$block_check{$block->{name}}=$block;
	      }
	  }
	@blocks = values %block_check;#create the non-redundant set
	#print out the blocks in order.  Note ones that are in reverse orientation
	foreach my $block (sort {$a->{match_start} <=> $b->{match_start} }@blocks)
	  {
	    next if $seen{$block->{name}};
	    push @$ordered2, {chr=>$block->{name}, rev=>$block->{rev}};
	    $seen{$block->{name}}=1;
	  }
      }
    return $ordered2;
#    ($ordered1, $ordered2) = ($ordered2, $ordered1) if $switched;
#    return $ordered1, $ordered2;
  }

sub process_syn_block
  {
    my $block = shift;
    my ($head, @block) = split/\n/, $block;
    my ($block_num, $score, $seq1, $seq2, $strand, $num_pairs) = 
      split/\t/, $head;
    my $rev = $head =~/r/ ? 1 : 0;
    my ($seq1_start, $seq1_stop, $seq2_start, $seq2_stop);
    #absolute start and stop can give rise to problems if the ends actually hit something far away from the rest of the sytnenic pairs.  Calculating the "mean" position will circumvent this problem
    my @start1;
    my @stop1;
    my @start2;
    my @stop2;
    foreach my $item (@block)
      {
	chomp $item;
	next unless $item;
	my @item = split /\t/, $item;
	next unless ($item[2] && $item[3] && $item[6] && $item[7]);
	push @start1, $item[2];
	push @stop1, $item[3];
	push @start2, $item[6];
	push @stop2, $item[7];
      }
    #remove the ends;
    @start1 = sort {$a<=>$b} @start1;
    @stop1 = sort {$a<=>$b} @stop1;
    @start2 = sort {$a<=>$b} @start2;
    @stop2 = sort {$a<=>$b} @stop2;
    shift @start1 if scalar(@start1) >3;
    pop @start1 if scalar(@start1) >2;
    shift @stop1 if scalar(@stop1) >3;
    pop @stop1  if scalar(@stop1) >2;
    shift @start2 if scalar(@start2) >3;
    pop @start2 if scalar(@start2) >2;
    shift @stop2 if scalar(@stop2) >3;
    pop @stop2  if scalar(@stop2) >2;
    map {$seq1_start+=$_} @start1;
    map {$seq1_stop+=$_} @stop1;
    map {$seq2_start+=$_} @start2;
    map {$seq2_stop+=$_} @stop2;
    $seq1_start = $seq1_start/scalar(@start1);
    $seq1_stop = $seq1_stop/scalar(@stop1);
    $seq2_start = $seq2_start/scalar(@start2);
    $seq2_stop = $seq2_stop/scalar(@stop2);
    my %seq1 = (
		name=>$seq1,
		start=>$seq1_start,
		stop=>$seq1_stop,
		match_start=>$seq2_start,
		match_stop=>$seq2_stop,
		score=>$score,
		rev=>$rev,
		);
    my %seq2 = (
		name=>$seq2,
		start=>$seq2_start,
		stop=>$seq2_stop,
		match_start=>$seq1_start,
		match_stop=>$seq1_stop,
		score=>$score,
		rev=>$rev,
		);
    return \%seq1, \%seq2;
  }
