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
    __PACKAGE__->mk_accessors('file', 'org1_ordered', 'org2_ordered',
    	'chrs1_scores', 'chrs2_scores', 'blocks1', 'blocks2', 'chrs1', 'chrs2',
    	'dsgid1', 'dsgid2');
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
	#s/#//g; 		# mdb removed 3/20/13 issue 55
	next if /^#/;	# mdb added 3/20/13 issue 55
        eval {
            my ($block1, $block2) = $self->process_syn_block($_);
            push @$blocks1 , $block1;
            push @$blocks2 , $block2;
        };

        warn "Unable to process syntenic block: $@" if $@;
      }
    close IN;
    $/="\n";

    #which organism has more pieces (aka chromosomes, contigs, etc.)  we are going to assemble the one with more pieces
#    print Dumper \@blocks;
    my $chrs1={};
    my $chrs2={};
    my $chrs1_scores={};
    my $chrs2_scores={};
    my $dsgid1;
    my $dsgid2;
    #need to assign a block to the highest scoring matching chromosomes
#    print  Dumper $blocks1, $blocks2;
    foreach my $item (@$blocks1)
#      foreach my $item (@$blocks1)
      {
	$chrs1->{$item->{name}}{$item->{match}}{count}++;
	$chrs1->{$item->{name}}{$item->{match}}{score}+=$item->{score};
	$chrs1->{$item->{name}}{$item->{match}}{hi_score}=$item->{score} unless $chrs1->{$item->{name}}{$item->{match}}{hi_score};
	$chrs1->{$item->{name}}{$item->{match}}{rev}+=$item->{score} if $item->{rev}; #needs to be a weighted metric based on this size of the block.  Bigger blocks in rev orientation get more weight
	$chrs1->{$item->{name}}{$item->{match}}{start}=$item->{match_start} unless $chrs1->{$item->{name}}{$item->{match}}{start};
	#update start position if new start is less than existing start AND the block higher scoring
	$chrs1->{$item->{name}}{$item->{match}}{start}=$item->{match_start} if
	    $chrs1->{$item->{name}}{$item->{match}}{hi_score} < $item->{score} ;
	$chrs1->{$item->{name}}{$item->{match}}{hi_score}=$item->{score} if $chrs1->{$item->{name}}{$item->{match}}{hi_score}<$item->{score};
	$chrs1_scores->{$item->{name}}+=$item->{score};
	$dsgid1 = $item->{dsgid};
      }
    foreach my $item (sort {$a->{name} cmp $b->{name} }@$blocks2)
      {
	my $print =0;
	$chrs2->{$item->{name}}{$item->{match}}{count}++;
	$chrs2->{$item->{name}}{$item->{match}}{score}+=$item->{score};
	$chrs2->{$item->{name}}{$item->{match}}{hi_score}=$item->{score} unless $chrs2->{$item->{name}}{$item->{match}}{hi_score};
	$chrs2->{$item->{name}}{$item->{match}}{rev}+=$item->{score} if $item->{rev}; #needs to be a weighted metric based on this size of the block.  Bigger blocks in rev orientation get more weight
	$chrs2->{$item->{name}}{$item->{match}}{start}=$item->{match_start} unless $chrs2->{$item->{name}}{$item->{match}}{start};
	#update start position if new start is less than existing start AND the block higher scoring
	$chrs2->{$item->{name}}{$item->{match}}{start}=$item->{match_start} if
	    $chrs2->{$item->{name}}{$item->{match}}{hi_score} < $item->{score} ;
	$chrs2->{$item->{name}}{$item->{match}}{hi_score}=$item->{score} if $chrs2->{$item->{name}}{$item->{match}}{hi_score}<$item->{score};
	$chrs2_scores->{$item->{name}}+=$item->{score};
	$dsgid2 = $item->{dsgid};
      }
    #blocks1 will contain fewer chromosomes; blocks2 will be ordered by it.
    my $switched =keys %$chrs1 > keys %$chrs2 ? 1 : 0;
    if ($switched)
      {
	($blocks1, $blocks2) = ($blocks2, $blocks1);
	($chrs1, $chrs2) = ($chrs2, $chrs1);
	($chrs1_scores, $chrs2_scores) = ($chrs2_scores, $chrs1_scores);
	($dsgid1,$dsgid2)= ($dsgid2, $dsgid1);
      }

    my %best_chr2contigs;
    foreach my $contig (keys %$chrs2)
      {
	my ($best_chr) = sort { $chrs2->{$contig}{$b}{score} <=> $chrs2->{$contig}{$a}{score} } keys %{$chrs2->{$contig}};
	push @{$best_chr2contigs{$best_chr}}, {%{$chrs2->{$contig}{$best_chr}}, name=>$contig};
      }
    my $ordered1 =[]; #storage for ordered chromosomes
    my $ordered2 =[]; #storage for ordered chromosomes
    #sort blocks for chr1 so that the highest scoring ones are first
    foreach my $chr1 (sort{$chrs1_scores->{$b} <=> $chrs1_scores->{$a}} keys %$chrs1_scores)
      {
	next unless $best_chr2contigs{$chr1};
	my %data;
	$data{chr}=$chr1;
	$data{matching_chr}=[map {$_->{name}} sort {$a->{start}<=>$b->{start}} @{$best_chr2contigs{$chr1}}];
	$data{rev}=0;
	$data{dsgid}=$dsgid1;
	#print out the blocks in order.  Note ones that are in reverse orientation
	foreach my $block (sort {$a->{start}<=>$b->{start}} @{$best_chr2contigs{$chr1}})
	  {
	    my $rev = 1 if $block->{rev} && $block->{rev}/$block->{score} >= 0.5;
	    push @$ordered2, {chr=>$block->{name}, rev=>$rev, matching_chr=>[$chr1], dsgid=>$dsgid2};
	  }
	push @$ordered1,\%data;
      }
    ($ordered1, $ordered2) = ($ordered2, $ordered1) if $switched;
    $self->org1_ordered($ordered1);
    $self->org2_ordered($ordered2);
    $self->chrs1($chrs1);
    $self->chrs2($chrs2);
    $self->chrs1_scores($chrs1_scores);
    $self->chrs2_scores($chrs2_scores);
    $self->blocks1($blocks1);
    $self->blocks1($blocks2);
    $self->dsgid1($dsgid1);
    $self->dsgid2($dsgid2);
    return $ordered1, $ordered2, $dsgid1, $dsgid2;
  }

sub process_syn_block
  {
    my $self = shift;
    my $block = shift;
    my ($head, @block) = split/\n/, $block;
    my ($block_num, $score, $seq1, $seq2, $strand) = split/\t/, $head;
    unless ($score)
      {
	$score = scalar @block;
	my @line1 = split/\t/, $block[0];
	my @line2 = split/\t/, $block[-1];

	$seq1 = $line1[0];
	$seq2 = $line1[4];
	$strand = $line1[6] < $line2[6] ? 0 : "r"; #start position will be smaller for last entry than first entry if line is reversed
      }
    my $rev = $strand =~/r/ ? 1 : 0;

    my ($seq1_start, $seq1_stop, $seq2_start, $seq2_stop);
    #absolute start and stop can give rise to problems if the ends actually hit something far away from the rest of the sytnenic pairs.  Calculating the "mean" position will circumvent this problem
    my @start1;
    my @stop1;
    my @start2;
    my @stop2;
    my $identity=0; #arbitrary place to add a numeric that sums the percent identities.  This value may be used to try to break a tie if a block matches two places with an equal score
    my $dsgid1;
    my $dsgid2;
    foreach my $item (@block)
      {
	chomp $item;
	next unless $item;
	next if $item =~ /^#/;
	my @item = split /\t/, $item;
	push @start1, $item[2] if $item[2] =~ /^\d+$/;
	push @stop1, $item[3] if $item[3] =~ /^\d+$/;
	push @start2, $item[6] if $item[6] =~ /^\d+$/;
	push @stop2, $item[7] if $item[7] =~ /^\d+$/;
	my @match = split /\|\|/,$item[1];
	$identity+=$match[8] if $match[8];
	unless ($dsgid1)
	  {
	    ($dsgid1) = split /_/, $item[0];
	    $dsgid1 =~ s/^a//;
	    $dsgid1 =~ s/^b//;
	    ($dsgid2) = split /_/, $item[4];
	    $dsgid2 =~ s/^a//;
	    $dsgid2 =~ s/^b//;
	  }
      }
    my $num_pairs = scalar @start1;
    #remove the ends;
    @start1 = sort {$a<=>$b} @start1;
    @stop1 = sort {$a<=>$b} @stop1;
    @start2 = sort {$a<=>$b} @start2;
    @stop2 = sort {$a<=>$b} @stop2;
    my $absolute_start1 = $start1[0];
    my $absolute_stop1 = $stop1[-1];
    my $absolute_start2 = $start2[0];
    my $absolute_stop2 = $stop2[-1];
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
		absolute_start=>$absolute_start1,
		stop=>$seq1_stop,
		absolute_stop=>$absolute_stop1,
		match_start=>$seq2_start,
		match_stop=>$seq2_stop,
		match_absolute_start=>$absolute_start2,
		match_absolute_stop=>$absolute_stop2,
		score=>$score,
		rev=>$rev,
		match=>$seq2,
		num_pairs=>$num_pairs,
		identity_sum=>$identity,
		dsgid=>$dsgid1,
		);
    my %seq2 = (
		name=>$seq2,
		start=>$seq2_start,
		stop=>$seq2_stop,
		absolute_start=>$absolute_start2,
		absolute_stop=>$absolute_stop2,
		match_start=>$seq1_start,
		match_stop=>$seq1_stop,
		match_absolute_start=>$absolute_start2,
		match_absolute_stop=>$absolute_stop1,
		score=>$score,
		rev=>$rev,
		match=>$seq1,
		num_pairs=>$num_pairs,
		identity_sum=>$identity,
		dsgid=>$dsgid2
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
