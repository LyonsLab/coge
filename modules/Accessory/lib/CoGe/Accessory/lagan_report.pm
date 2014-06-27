package CoGe::Accessory::lagan_report;

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
__PACKAGE__->mk_accessors('file', 'hsps', 'hsp_count', 'max_gap', 'query',
	'subject', 'qlength', 'slength', 'length_cutoff', 'percent_cutoff' );

###############################################################################
# lagan_report -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->hsp_count(0) unless $self->hsp_count;
	$self->length_cutoff(0) unless $self->length_cutoff;
	$self->percent_cutoff(0) unless $self->percent_cutoff;
	$self->max_gap(0) unless $self->max_gap;

	# init
	$self->process_file();
	return $self;
      }

sub process_file
  {

    my $self = shift;
    my $file = shift || $self->file;
    my @names;
    $/ = "\n>";
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
    {
     $_ =~ s/\>//g;
     my ($n, $seq) = split /\n/, $_, 2;
     $seq =~ s/\n//g;
     push @names,{name=>$n, seq=>$seq};
    }
    close IN;
    $/ = "\n";
    $self->query($names[0]{name});
    $self->subject($names[1]{name});
    my $ql = $names[0]{seq} =~ tr/ATGC/ATGC/;
    my $sl = $names[1]{seq} =~ tr/ATGC/ATGC/;
    $self->qlength($ql);
    $self->slength($sl);
    $self->_parseReport($names[0]{seq},$names[1]{seq});
    return $self;
}

sub _parseReport {
	my $self = shift;
	my $sq1 = shift;
	my $sq2 = shift;
	my @seq1 = make_array($sq1);
	my @seq2 = make_array($sq2);
	my @matches = make_align(\@seq1,\@seq2);
	my @hsps;
	my $hsp_count = $self->hsp_count();
	my $check=0;
	my ($gap1,$gap2) = (0,0);
	my ($stop1,$stop2)=(1,1);
	my ($align1,$align2,$align);
	my ($st1,$st2,$tmp1,$tmp2,$break,$gap_length)=(0,0,0,0,0,0);
	for (my $i=0;$i<scalar (@seq1);$i++)
	{
  	 $align1 .= $seq1[$i];
 	 $align2 .= $seq2[$i];
 	 $align .= $matches[$i];
  	 if ($seq1[$i] =~/-/)
   	  {$gap1++;$stop1--;$tmp2++;}
 	 if ($seq2[$i] =~/-/)
   	  {$gap2++;$stop2--;$tmp1++;}
   #print "$stop1\n";
  	 if ($gap1 && $seq1[$i] !~/-/)
  	  {$check = 1;}
  	 if ($gap2 && $seq2[$i] !~/-/)
  	  {$check = 2;}
   #print "check: $check\n";
  	 if ($check)
  	 {
   	  if ($check == 1)
    	  {
     	   if ($gap1 > $self->max_gap)
      	    {$gap_length = $gap1+1;}
      	   else
       	   {
            $gap1 = 0;
            $tmp2 = 0;
            $break = 1;
       	   }
          }
    	 else
     	 {
      	  if ($gap2 > $self->max_gap)
       	   {$gap_length = $gap2+1;}
       	  else
       	  {
           $gap2 = 0;
           $tmp2 = 0;
           $break = 1;
       	  }
     	 }

     	unless ($break)
        {
    	 $align1 = substr($align1,0,(length $align1) - $gap_length);
     	 $align2 = substr($align2,0,(length $align2) - $gap_length);
   	 $align = substr($align,0,(length $align) - $gap_length);
    	 if ($align =~/\|+/)
    	 {
          #if (($length >$min_align_length)&&($ident1 >=$min_ident)&&($ident2 >= $min_ident))
           my $hsp = $self->_processHSP($align1,$align2,$align,$tmp1,$stop1,$tmp2,$stop2);

           if ($hsp)
            {
	      $hsp_count++;
	      push @hsps, $hsp;
	      $self->hsp_count($hsp_count);
	    }
	}
          ($gap1,$gap2) = (0,0);
         $align1 = $seq1[$i];
         $align2 = $seq2[$i];
         $align = $matches[$i];
        }
        ($break,$check,$tmp1,$tmp2) = (0,0,0,0);
       }
     $stop1++;$stop2++;
    }
	$align1 = substr($align1,0,length $align1);
	$align2 = substr($align2,0,length $align2);
	$align = substr($align,0,length $align);
	if ($align =~/\|+/)
	  {
	    #if (($length >$min_align_length)&&($ident1 >=$min_ident)&&($ident2 >= $min_ident))
	    my $hsp = $self->_processHSP($align1,$align2,$align,$tmp1,$stop1,$tmp2,$stop2);
	    if ($hsp)
	      {
		$hsp_count++;
		push @hsps, $hsp;
		$self->hsp_count($hsp_count);
	      }
	  }
    $self->hsps(\@hsps);
    return $self;
}

sub _processHSP {
	my $self = shift;
	my $align1 = shift;
	my $align2 = shift;
	my $align = shift;
	my $tmp1 = shift;
	my $stop1 = shift;
	my $tmp2 = shift;
	my $stop2 = shift;
	my $align_length = length $align;
      	my ($ident1,$qm,$length1,$gaps1) = check_num_of_aligns($align1,$align);
      	my ($ident2,$sm,$length2,$gaps2) = check_num_of_aligns($align2,$align);
	$stop1=(($stop1-1)-$tmp1);my $start1 = (($stop1-$length1)+1);
	$stop2=(($stop2-1)-$tmp2);my $start2 = (($stop2-$length2)+1);
	while ($align !~ /^\|/)
    	   {$align1 =~s/^.//;$align2 =~s/^.//;$align =~s/^.//;$start1++;$start2++;}
    	while ($align !~ /\|$/)
       	   {$align1 =~s/.$//;$align2 =~s/.$//;$align =~s/.$//;$stop1--;$stop2--;}
	my $hsp_count = $self->hsp_count();
	$hsp_count++;

	my $hsp = new CoGe::Accessory::parse_report::HSP
           ({
    	     qpercent_id=>$ident1,
    	     spercent_id=>$ident2,
    	     match=>$qm,
	     percent_id=>$ident1,
    	     length=>$align_length,
    	     query_start=>$start1,
    	     query_stop=>$stop1,
    	     subject_start=>$start2,
    	     subject_stop=>$stop2,
    	     query_alignment=>$align1,
    	     subject_alignment=>$align2,
    	     alignment=>$align,
    	     query_gaps=>$gaps1,
    	     subject_gaps=>$gaps2,
    	     number=>$hsp_count,
	     pval=>"N/A",
	     strand=>1,
            });
        return 0 if $align_length < $self->length_cutoff || $self->percent_cutoff > $ident1 || $self->percent_cutoff > $ident2;
	return $hsp;
}

sub make_array
{
 my $str = shift;
 my @seq;
 foreach (split //,$str)
  {
   push @seq,$_;
  }
  return @seq;
} my ($gap1,$gap2) = (0,0);

 sub make_align
{
  my ($ref1,$ref2) = @_;
  my $tmp;
  my @align;
  for (my $i=0;$i<scalar(@$ref1);$i++)
  {
    $tmp = (($ref1->[$i] eq $ref2->[$i])&&($ref1->[$i] !~/-/)) ? "|" : " ";
    push @align,$tmp;
   }
  return @align;
}

sub check_num_of_aligns
{
  my $seq = shift;
  my $align = shift;
  my ($count1,$count2,$gaps) = (0,0,0);
  $count1 = $seq =~tr/ATGC/ATGC/;
  #print "\$count1 is $count1\n";
  $gaps = (length $seq) - $count1;
  $count2 = $align =~tr/\|/\|/;
  my $percent = sprintf("%.4f", ($count2/$count1)) * 100;
  #$gaps = 0 unless $gaps;
  return $percent,$count2,$count1,$gaps;
}

1;

__END__
