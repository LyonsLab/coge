package CoGe::Accessory::dialign_report;

use strict;
use warnings;
use CoGe::Accessory::parse_report::HSP;
use base qw(Class::Accessor);

use Data::Dumper;

BEGIN
  {
    use vars qw($VERSION);
    $VERSION = "0.01";
  }
__PACKAGE__->mk_accessors('file', 'hsps', 'hsp_count', 'query', 'subject',
	'qlength', 'slength', 'min_score', 'max_gap', 'split_score');

###############################################################################
# dialign_report -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->split_score(0) unless $self->split_score;
	$self->hsp_count(0) unless $self->hsp_count;
	$self->min_score(1) unless $self->min_score;
	$self->max_gap(2) unless $self->max_gap;
	$self->process_file();
	return $self;
      }

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    my @data;
    $/ = "=\n";
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
    {
     push @data,$_;
    }
    close IN;
    $/ = "\n";
    splice(@data,0,1);
    $self->parse_header($data[0]);
    $self->parse_seq($data[1]);
    return $self;
}

sub parse_header
{
 my $self = shift;
 my $str = shift;
 my ($clean_str,$nm1,$ln1,$nm2,$ln2) = ("",0,0,0,0);
 foreach my $bits (split /\n/,$str)
 {
   $clean_str .= "$bits" unless $bits !~ /\d+/;
  }
  if ($clean_str =~/^\s*\d.\s(\w+\.?\d*)\s+(\d+)\s+\d.\s(\w+\.?\d*)\s+(\d+)/)
  {($nm1,$ln1,$nm2,$ln2) = ($1,$2,$3,$4);}
  $self->query($nm1);
  $self->qlength($ln1);
  $self->subject($nm2);
  $self->slength($ln2);
  return $self;
}

sub parse_seq
{
  my $self = shift;
  my $str = shift;
  my ($sq1,$sq2,$nums);
  my (@holder,@hsps);
  foreach my $parts (split /\n/,$str)
  {
   push @holder,$parts unless $parts !~ /\d+/;
  }
  for (my $i=0;$i<scalar(@holder);$i+=3)
  {
    $sq1 .= $holder[$i];
    $sq2 .= $holder[$i+1];
    $nums .= $holder[$i+2];
  }
  $nums =~ s/\D//g;
  my $nm1 = $self->query;
  my $nm2 = $self->subject;
  $sq1 =~ s/$nm1//g;
  $sq1 =~ s/(\s|\d)//g;
  $sq2 =~ s/$nm2//g;
  $sq2 =~ s/(\s|\d)//g;
  my ($scores) = $self->get_score_aligns($nums);
  #print Dumper @$scores;
  my ($seq1,$seq2) = $self->get_align($scores,$sq1,$sq2);
  #print Dumper @$seq1;
  #print Dumper @$seq2;
  for (my $i=0;$i<scalar(@$seq1);$i++)
  {
    my $hsp = $self->_getHSP(%{$seq1->[$i]},%{$seq2->[$i]},%{$scores->[$i]});
    push @hsps,$hsp if $hsp;
  }
  $self->hsps(\@hsps);
  return $self;
}

sub _getHSP
{
  my $self = shift;
  my %opts = @_;
  #print Dumper %opts;
  my $sq1 = $opts{align1};
  my $sq2 = $opts{align2};
  my $start1 = $opts{start1};
  my $start2 = $opts{start2};
  my $stop1 = $opts{stop1};
  my $stop2 = $opts{stop2};
  my $scores = $opts{score};
  #print "scores are $scores\n";
  #my $stop = $opts{stop};
  #my $start = $opts{start};
  my $align = $self->get_pipes($sq1,$sq2,$scores);
  #my $align = "|||";
  my $match = $align =~ tr/\|/\|/;
  while ($align && $align !~ /^\|/)
       {$sq1 =~s/^.//;$sq2 =~s/^.//;$align =~s/^.//;$start1++;$start2++;}
  while ($align && $align !~ /\|$/)
       {$sq1 =~s/.$//;$sq2 =~s/.$//;$align =~s/.$//;$stop1--;$stop2--;}
  my $hsp_count = $self->hsp_count;
  $hsp_count++;
  my $length = length $sq1;
  my $hsp = new CoGe::Accessory::parse_report::HSP
        ({
   	  score=>$scores,
   	  strand=>1,
    	  match=>$match,
    	  length=>$length,
    	  query_start=>$start1,
    	  query_stop=>$stop1,
    	  subject_start=>$start2,
    	  subject_stop=>$stop2,
    	  query_alignment=>$sq1,
    	  subject_alignment=>$sq2,
    	  alignment=>$align,
    	  number=>$hsp_count,
        });
    $self->hsp_count($hsp_count) if $hsp;
    return $hsp;
}

sub get_align
{
 my $self = shift;
 my ($vals) = shift;
 my $sq1 = shift;
 my $sq2 = shift;
 #print Dumper @$vals;
 my ($start,$stop,$score,$length);
 my ($start1,$stop1,$start2,$stop2);
 my ($tmp1,$tmp2);
 my $seq1 = [];
 my $seq2 = [];
 for (my $i=0;$i < scalar(@$vals);$i++)
 {
   $start = $vals->[$i]{start};
   #print "start is $start\n";
   $stop = $vals->[$i]{stop};
   #print "stop is $stop\n";
   $score = $vals->[$i]{score};
   $length = (length $score) - 1;
   $tmp1 = substr($sq1,$start,$length);
   #print "tmp1 is $tmp1\n";
   $tmp2 = substr($sq2,$start,$length);
   #print "tmp2 is $tmp2\n";
   ($start1,$stop1) = $self->counting(substr($sq1,0,$stop),$tmp1);
   ($start2,$stop2) = $self->counting(substr($sq2,0,$stop),$tmp2);
   push @$seq1,{start1=>$start1,stop1=>$stop1,align1=>$tmp1};
   push @$seq2,{start2=>$start2,stop2=>$stop2,align2=>$tmp2};
 }
 return $seq1,$seq2;
}

# sub get_align
# {
#   my $self = shift;
#   my $str = shift;
#   my ($count,$hit,$align,$start,$stop) = (0,0,"",0,0);
#   my (@seq,@aligns,@starts,@stops);
#   foreach my $char (split //,$str)
#   {
#    push @seq,$char if $char =~ /[ATGC-]/i;
#   }
#   for(my $i=0;$i<scalar(@seq);$i++)
#   {
#       $count++ if $seq[$i] !~ /-/;
#       if ($seq[$i] =~/[ATGC]/)
#       {
#         $align .= $seq[$i];
#         $hit++;
#         if ($seq[$i+1] !~ /[ATGC]/)
#         {
#          $stop = $count;
#          $start = ($stop - $hit)+1;
#          push @starts,$start;
#          push @stops,$stop;
#          push @aligns,$align;
#          ($start,$stop,$align,$hit) = (0,0,"",0);
#         }
#       }
#    }
#
#
#    return (\@starts,\@stops,\@aligns);
#}

sub get_pipes
{
  my $self = shift;
  my ($seq1,$seq2,$scores) = @_;
  #print "($seq1,$seq2,$scores)\n";
  my $align;
  my $length = length $seq1;
  while ($seq1)
  {
    $align .= ((substr($seq1,0,1) eq substr($seq2,0,1)) && substr($scores,0,1)) ? "|" : " ";
    $seq1 =~ s/^.//;
    $seq2 =~ s/^.//;
    $scores =~ s/^.//;
  }
  return $align;
}

sub get_score_aligns
{
  my $self = shift;
  my $num = shift;
  my $info = [];
  my $min_score = $self->min_score;
  my $max_gap = $self->max_gap;
  my $split_scores = $self->split_score;
  my $count = 0;
  my $check;
  my $start = 0;
  my $stop = 0;
  my $loc = 0;
  my $hits = 0;
  my $tmp = "";
  foreach my $digits (split //,$num)
  {
   next unless $digits =~ /\d/;
   #print $digits;
   $count++;
   unless ($tmp)
    {next if $digits == 0;}
   if ($split_scores)
   {
     $check = $1 if $tmp=~/(\d$)/;
     if ($tmp && ($check != $digits))
     {
       #print "regex if statement worked! $check is not $digits\n";
       #print "tmp is $tmp\n";
       $stop = $count - 1;
       $start = $stop - (length $tmp);
       #print "stop: $stop\n";
	#print "start: $start\n";
       push @$info,{start=>$start,stop=>$stop,score=>$tmp};
       $tmp = $digits unless $digits == 0;
       $tmp = "" if $digits == 0;
     }
     else
      {$tmp .= $digits;}
    }
    else
    {
   #print "$count\n";
	if ($digits < $min_score)
	{
	$hits++;
	#print "hits: $hits\n";
	if ($hits > $max_gap)
	{
	$hits--;
	#print "score before edit: $tmp\n";
	$tmp = substr($tmp,0,((length $tmp) - $hits));
	#print "score after edit: $tmp\n";
	$stop = $count - $hits;
	#print "stop: $stop\n";
	#print "length of tmp is",length $tmp,"\n";
	$start = $stop - (length $tmp);
	# print "start: $start\n";
	push @$info,{start=>$start,stop=>$stop,score=>$tmp};
	($tmp,$start,$hits) = ("",0,0);
	}
	else
	{$tmp .= $digits;}
	}
	else
	{
	$hits = 0 if $hits;
	$tmp .= $digits;
	}
    }
  }
  if ($tmp)
  {
    while ($tmp=~/(\d)$/)
    {
      if ($1 < $min_score)
      {
       $tmp=~s/.$//;
       $count--;
      }
      else
       {last;}
     }
    $start = ($count - (length $tmp));
    push @$info,{start=>$start,stop=>$count,score=>$tmp};
    }
  return $info;
}

sub counting
{
  my $self = shift;
  my $seq = shift;
  my $align = shift;
  $align =~ s/-//g;
  my $length = length $align;
  my $count = 0;
  foreach my $char (split //,$seq)
  {
    $count++ if $char =~ /[atgc]/i;
  }
  my $start = ($count - $length)+1;
  return $start,$count;
}

1;

__END__
