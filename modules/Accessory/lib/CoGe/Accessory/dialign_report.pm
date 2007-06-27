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
__PACKAGE__->mk_accessors qw(file hsps hsp_count query subject qlength slength);

###############################################################################
# dialign_report -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->hsp_count(0) unless $self->hsp_count;;
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
  $nums =~ s/\s//g;
  $sq1 =~ s/$self->query//g;
  $sq2 =~ s/$self->subject//g;
  my ($start1,$stop1,$align1) = $self->get_align($sq1);
  my ($start2,$stop2,$align2) = $self->get_align($sq2);
  my ($values) = $self->get_values($nums);
  for (my $i=0;$i<scalar(@$start2);$i++)
  {
    my $hsp = $self->_getHSP($start1->[$i],$stop1->[$i],$align1->[$i],$start2->[$i],$stop2->[$i],$align2->[$i],$values->[$i]);
    push @hsps,$hsp if $hsp;
  }
  $self->hsps(\@hsps);
  return $self;
}

sub _getHSP
{
  my $self = shift;
  my ($start1,$stop1,$sq1,$start2,$stop2,$sq2,$values) = @_;
  my $align = $self->get_pipes($sq1,$sq2);
  my $match = $align =~ tr/|/|/;
  while ($align !~ /^\|/)
        {$sq1 =~s/^.//;$sq2 =~s/^.//;$align =~s/^.//;}
  while ($align !~ /\|$/)
       	{$sq1 =~s/.$//;$sq2 =~s/.$//;$align =~s/.$//;}
  my $hsp_count = $self->hsp_count;
  $hsp_count++;
  my $length = length $sq1;
  my $hsp = new CoGe::Accessory::parse_report::HSP
        ({
   	  score=>$values,
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
  my $str = shift;
  my ($hit,$align,$start,$stop) = (0,"",0,0);
  my (@seq,@aligns,@starts,@stops);
  foreach my $char (split //,$str)
  {
   push @seq,$char if $char =~ /[ATGC]/i;
  }
  for(my $i=0;$i<scalar(@seq);$i++)
  {
      if ($seq[$i] =~/[ATGC]/)
      {
        $align .= $seq[$i];
        $hit++;
        if ($seq[$i+1] !~ /[ATGC]/)
        {
         $stop = $i+1;
         $start = $stop - $hit;
         push @starts,$start;
         push @stops,$stop;
         push @aligns,$align;
         ($start,$stop,$align) = (0,0,"");
        }
       }
   }
   return (\@starts,\@stops,\@aligns);
}

sub get_pipes
{
  my $self = shift;
  my ($seq1,$seq2) = @_;
  my $align;
  for (my $i=0;$i<(length $seq1);$i++)
  {
    $align .= $seq1=~/^./ eq $seq2=~/^./ ? "|" : " ";
    $seq1 =~ s/^.//;
    $seq2 =~ s/^.//;
  }
  return $align;
}
 
sub get_values
{
  my $self = shift;
  my $num = shift;
  my @nums;
  foreach my $digits (split /0+/,$num)
  {
   push @nums,$digits if $digits;
  }
  return \@nums;
}

1;

__END__
