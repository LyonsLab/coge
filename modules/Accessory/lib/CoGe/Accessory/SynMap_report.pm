package CoGe::Accessory::SynMap_report;
use strict;
use base qw(Class::Accessor) ;
use Data::Dumper ;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $AUTOLOAD);
    $VERSION     = 0.1;
    @ISA         = qw (Exporter Class::Accessor );
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->mk_accessors(qw(file org1_ordered org2_ordered)) ;
}


#this sub is used to try to get a build order of contigs for WGS data against a reference genome
#given a set of scaffolds, it will order them such that make an ordered syntenic path along the reference genome
sub parse_syn_blocks
  {
    my $self = shift;
    my %opts = @_;
    my $file = $opts{file};
    $file = $self->file() unless $file;
    $self->file($file);
    my $blocks1=[];
    my $blocks2=[];
    open (IN, $file) || die "Can't open $file for reading: $!";
    $/ = "\n#";
    while (<IN>) #get blocks
      {
	next unless $_;
	s/#//g;
	my ($block1, $block2) = $self->process_syn_block($_);
	push @$blocks1 , $block1;
	push @$blocks2 , $block2;
      }
    close IN;
    $/="\n";

    #which organism has more pieces (aka chromosomes, contigs, etc.)  we are going to assemble the one with more pieces
#    print Dumper \@blocks;
    my $chrs1={};
    my $chrs2={};
    my $chrs1_scores={};
    my $chrs2_scores={};
    #need to assign a block to the highest scoring matching chromosomes

    foreach my $item (sort {$a->{name} cmp $b->{name} || $a->{match} cmp $b->{match} || $a->{score} <=> $b->{score} || $a->{identity_sum} <=> $b->{identity_sum}} @$blocks1)
#      foreach my $item (@$blocks1)
      {
	$chrs1->{$item->{name}}{$item->{match}}{count}++;
	$chrs1->{$item->{name}}{$item->{match}}{score}+=$item->{score};
	$chrs1_scores->{$item->{name}}+=$item->{score};
      }
    foreach my $item (@$blocks2)
      {
	$chrs2->{$item->{name}}{$item->{match}}{count}++;
	$chrs2->{$item->{name}}{$item->{match}}{score}+=$item->{score};
	$chrs2_scores->{$item->{name}}+=$item->{score};
      }
    #blocks1 will contain fewer chromosomes; blocks2 will be ordered by it.
    my $switched =keys %$chrs1 > keys %$chrs2 ? 1 : 0;
    if ($switched)
      {
	($blocks1, $blocks2) = ($blocks2, $blocks1);
	($chrs1, $chrs2) = ($chrs2, $chrs1);
	($chrs1_scores, $chrs2_scores) = ($chrs2_scores, $chrs1_scores);
      }


    #need to use only the highest scoring blocks in set2 to take into account duplications
    my %best_blocks2;
    foreach my $block (@$blocks2)
      {
	# is this chromosome already assigned to another match?  If so, is the one we have better in terms of having a higher score or (same score and higher identity_sum)?  If so, let's use that one.
	if ($best_blocks2{$block->{name}})
	  { 
	    if ( $block->{score} > $best_blocks2{$block->{name}}->{score} ||
		( $block->{score} == $best_blocks2{$block->{name}}->{score} && $block->{identity_sum} > $best_blocks2{$block->{name}}->{identity_sum})
	       )
	      {
		$best_blocks2{$block->{name}} = $block;
	      }
	  }
	else
	  {
	    $best_blocks2{$block->{name}} = $block;
	  }
      }
    $blocks2 = [values %best_blocks2];


    my $ordered1 =[]; #storage for ordered chromosomes
    my $ordered2 =[]; #storage for ordered chromosomes
    #sort blocks for chr1 so that the highest scoring ones are first
    foreach my $chr1 (sort{$chrs1_scores->{$b} <=> $chrs1_scores->{$a}} keys %$chrs1_scores)
      {
	push @$ordered1, {chr=>$chr1};
	my @blocks;
	foreach my $block (@$blocks2)
	  {
	    push @blocks, $block if $block->{match} eq $chr1;
	  }
	#Need to check if a given chromosome in @blocks occurs more than once.  This will happen if there is a segmental duplication, or something along those lines.
	my %block_check;
	foreach my $block (@blocks)
	  {
	    if ($block_check{$block->{name}})
	      {
		$block_check{$block->{name}}=$block if $block->{num_pairs} > $block_check{$block->{name}}{num_pairs}; #more pairs
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
	    push @$ordered2, {chr=>$block->{name}, rev=>$block->{rev}, match=>$block->{match}};
	  }
      }
    ($ordered1, $ordered2) = ($ordered2, $ordered1) if $switched;
    $self->org1_ordered($ordered1);
    $self->org2_ordered($ordered2);
    return $ordered1, $ordered2;
  }


sub process_syn_block
  {
    my $self = shift;
    my $block = shift;
    my ($head, @block) = split/\n/, $block;
    my ($block_num, $score, $seq1, $seq2, $strand) = 
      split/\t/, $head;
    my $rev = $strand =~/r/ ? 1 : 0;

    my ($seq1_start, $seq1_stop, $seq2_start, $seq2_stop);
    #absolute start and stop can give rise to problems if the ends actually hit something far away from the rest of the sytnenic pairs.  Calculating the "mean" position will circumvent this problem
    my @start1;
    my @stop1;
    my @start2;
    my @stop2;
    my $identity=0; #arbitrary place to add a numeric that sums the percent identities.  This value may be used to try to break a tie if a block matches two places with an equal score
    foreach my $item (@block)
      {
	chomp $item;
	next unless $item;
	my @item = split /\t/, $item;
	push @start1, $item[2];
	push @stop1, $item[3];
	push @start2, $item[6];
	push @stop2, $item[7];
	my @match = split /\|\|/,$item[1];
	$identity+=$match[8] if $match[8];
      }
    my $num_pairs = scalar @start1;
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
    $seq1 =~ s/.*?_//;
    $seq2 =~ s/.*?_//;
    my %seq1 = (
		name=>$seq1,
		start=>$seq1_start,
		stop=>$seq1_stop,
		match_start=>$seq2_start,
		match_stop=>$seq2_stop,
		score=>$score,
		rev=>$rev,
		match=>$seq2,
		num_pairs=>$num_pairs,
		identity_sum=>$identity,
		);
    my %seq2 = (
		name=>$seq2,
		start=>$seq2_start,
		stop=>$seq2_stop,
		match_start=>$seq1_start,
		match_stop=>$seq1_stop,
		score=>$score,
		rev=>$rev,
		match=>$seq1,
		num_pairs=>$num_pairs,
		identity_sum=>$identity,
		);
    return \%seq1, \%seq2;
  }



1;

########################################### Main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

SynMap_report

=head1 SYNOPSIS

  use CoGe::Accessory::SynMap_report



=head1 DESCRIPTION

This module is for parsing the output of SynMap

=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

Eric Lyons
iPlant/UC Berkeley
elyons.berkeley@gmail.com

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

############################################# main pod documentation end ##

