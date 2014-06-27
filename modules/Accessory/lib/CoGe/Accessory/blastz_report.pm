package CoGe::Accessory::blastz_report;

use strict;
use POSIX;
use CoGe::Accessory::blastz_report::seq;
use CoGe::Accessory::blastz_report::segment;
use CoGe::Accessory::parse_report::HSP;
use base qw(Class::Accessor);
use Carp;
use Data::Dumper;

BEGIN
  {
    use vars qw($VERSION);
    $VERSION = "0.01";
  }
__PACKAGE__->mk_accessors('hsps', 'hsp_count', 'file', '_params', 'seq_info',
	'file_base', 'query', 'subject', 'match');

###############################################################################
# blastz report parser -- Eric Lyons All rights and lefts reserved!
###############################################################################
#################### subroutine header begin ####################

=head2 new

 Usage     : my $blastz_report = new CoGe::Accessory::blastz_report($blastz_file);
 Purpose   : parse a blastz report, generate this object
 Returns   : a blastz_report object
 Argument  : a string that points to a blastz file
 Throws    : lots of errors and warnings if the blastz file is missing or can't be read
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub new
  {
    my $proto = shift;
    my $opts = shift;
    $opts = {} unless $opts;
    my $class = ref($proto) || $proto;
    my $self = bless ({%$opts}, $class);
    $self->hsp_count(0);
    $self->hsps([]) unless $self->hsps;
    $self->process_file();
    unless ($self->file_base())
      {
	my ($base) = $self->file =~ /^(.*?)[^\/]*$/;
	$base = "./".$base unless $base =~ /^\//;
	$self->file_base($base);
      }
    return $self;
}

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    unless (-r $file)
      {
	warn "can't read $file in blastz_report.pm.  nothing to process\n";
	return;
      }
    open (IN, $file);
    $/ = "}\n";
    while (<IN>)
      {
	s/#.*\n//g;
	$self->parse_params($_) if /^d/;
	$self->parse_file_names($_) if /^s/;
	$self->parse_fasta_header($_) if /^h/;
	$self->parse_alignment($_) if /^a/;
      }
    $self->parse_sequences(); #let's see if we get lucky with finding sequence files
    close IN;
  }

sub parse_params
  {
    my $self = shift;
    my $data = shift;
    my %results;
    my $matrix = 0;
    my %matrix = (1=>"A",
		  2=>"C",
		  3=>"G",
		  4=>"T",
		  );
    foreach my $line (split /\n/, $data)
      {
	next if $line =~ /{|}/; #skip brackets
	if ($line =~ /=/)
	  {
	    $matrix = 0;
	    foreach (split/,/, $line)
	      {
		s/\s+|"//g;
		my ($param, $val) = split /=/;
		$results{$param}=$val;
	      }
	  }
	if ($line =~ /\s+A\s+C\s+/i) #ATCG id line
	  {
	    $matrix = 1;
	    next;
	  }
	if ($matrix)
	  {
	    $line =~ s/^\s+//;
	    my @line = split /\s+/, $line;
	    for (1..4)
	      {
		$results{matrix}{$matrix{$matrix}}{$matrix{$_}} = $line[$_-1];
	      }
	    $matrix++;
	  }

      }
    $self->_params(\%results);
  }
#################### subroutine header begin ####################

=head2 matrix string

 Usage     : my $matrix = $blastz->matrix_string();
 Purpose   : generates a string of the matrix used
 Returns   : a string
 Argument  : none (blastz file is parsed during object creation)
 Throws    :
 Comment   : Example file:
        A       C       G       T
A       91      -114    -31     -123
C       -114    100     -125    -31
G       -31     -125    100     -114
T       -123    -31     -114    91

See Also   :

=cut

#################### subroutine header end ####################

sub matrix_string
  {
    my $self=shift;
    my $results;
    if (my $matrix = $self->matrix)
      {
	$results = "\t".join("\t", sort keys %$matrix)."\n";
	foreach my $i (sort keys %$matrix)
	  {
	    $results .= join ("\t", $i, map {$matrix->{$i}{$_}} sort keys %{$matrix->{$i}})."\n";
	  }
      }
    return $results;
  }

#################### subroutine header begin ####################

=head2 matrix

 Usage     : my $matrix = $blastz->matrix();
 Purpose   : gets a hash of the matrix used in blastz
 Returns   : a hashref
 Argument  : none
 Throws    : carps if matrix has not been parsed
 Comment   : hash structure:
$VAR1 = {
          'A' => {
                   'A' => '91',
                   'T' => '-123',
                   'C' => '-114',
                   'G' => '-31'
                 },
          'T' => {
                   'A' => '-123',
                   'T' => '91',
                   'C' => '-31',
                   'G' => '-114'
                 },
          'C' => {
                   'A' => '-114',
                   'T' => '-31',
                   'C' => '100',
                   'G' => '-125'
                 },
          'G' => {
                   'A' => '-31',
                   'T' => '-114',
                   'C' => '-125',
                   'G' => '100'
                 }
        };

See Also   :

=cut

#################### subroutine header end ####################
sub matrix
  {
    my $self = shift;
    if (ref ($self->_params) =~ /hash/i && $self->_params->{matrix})
      {
	return $self->_params->{matrix};
      }
    else
      {
	my $file = $self->file ||"No file!";
	carp "Ain't got no matrix!  Failed blastz parsing perhaps?  Missing matrix?  Check file: ".$file."\n";
      }
  }
#################### subroutine header begin ####################

=head2 params

 Usage     : my $params = $self->params();
 Purpose   : get a hashref of the parameters used in the blastz run
 Returns   : a hashref
 Argument  : none
 Throws    : carps if no parameters exist
 Comment   : Example hashref:

$VAR1 = {
          'O' => '400',
          'M' => '50',
          'K' => '3000',
          'L' => '3000',
          'E' => '30'
        };

See Also   :

=cut

#################### subroutine header end ####################

sub params
  {
    my $self = shift;
    my %params;
    if (ref ($self->_params) =~ /hash/i)
      {
	foreach my $item (sort keys %{$self->_params})
	  {
	    next if $item =~ /matrix/;
	    $params{$item} = $self->_params->{$item};
	  }
      }
    else
      {
	my $file = $self->file ||"No file!";
	carp "Ain't got no params!  Failed blastz parsing perhaps?  Missing matrix?  Check file: ".$file."\n";
      }
    return \%params;
  }
#################### subroutine header begin ####################

=head2 param_string

 Usage     : my $pstring = $blastz->param_string();
 Purpose   : gets a string of the blastz parameters (without matrix)
 Returns   : a string
 Argument  : none
 Throws    : carps through ->params() if there are no parameters
 Comment   : example string
           : E=30, K=3000, L=3000, M=50, O=400

See Also   :

=cut

#################### subroutine header end ####################

sub param_string
  {
    my $self = shift;
    my $params = $self->params;
    return join (", ", map {$_."=". $params->{$_}} sort keys %$params);
  }

sub parse_file_names
  {
    my $self = shift;
    my $data = shift;
    my %results;
    my $seq_count = 1;
    foreach my $line (split /\n/, $data)
      {
	next if $line =~ /{|}/; #skip brackets
	$line =~ s/"//g;
	next unless $line;
	$line =~ s/^\s+//;
	my ($name, $start, $stop, $strand, $contig) = split /\s+/, $line;
	$name =~ s/-$//; #remove "-" from end of file name if blastz auto-used reverse complement
	$strand = $strand ? "-1" : "1";
	$results{$seq_count} = new CoGe::Accessory::blastz_report::seq({
								       file=>$name,
								       start=>$start,
								       stop=>$stop,
								       length=>($stop-$start+1),
								       strand=>$strand,
								       contig=>$contig,
								      });
	$seq_count++;
      }
    $self->seq_info(\%results);
  }

sub parse_fasta_header
  {
    my $self = shift;
    my $data = shift;
    my $seq_count = 1;
    foreach my $line (split /\n/, $data)
      {
	next if $line =~ /{|}/; #skip brackets
	$line =~ s/"//g;
	next unless $line;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//g;
	$line =~ s/^>//;
	$self->seq_info->{$seq_count}{header}=$line;
	$line =~ s/>//;
	if ($self->query)
	  {
	    $self->subject($line);
	  }
	else
	  {
	    $self->query($line);
	  }
	$seq_count++;
      }
  }
#################### subroutine header begin ####################

=head2 seq1_info

 Usage     : my $seq = $blastz->seq1_info();
 Purpose   : get sequence information for seq1 from a blastz run
 Returns   : a CoGe::Accessory::blastz_report::seq object
 Argument  : none
 Throws    :
 Comment   : the seq object has accessor functions for:
              file => file name
              start => start position
              stop => stop position
              length => length of sequence
              strand => strand used
              contig => something blastz uses when comparing contigs.  Not sure what this means, check blastz docs.

See Also   :

=cut

#################### subroutine header end ####################

sub seq1_info
  {
    my $self = shift;
    if (ref ($self->seq_info) =~ /hash/i)
      {
	return ($self->seq_info->{1});
      }
  }
#################### subroutine header begin ####################

=head2 seq2_info

 Usage     : my $seq = $blastz->seq2_info();
 Purpose   : get sequence information for seq2 from a blastz run
 Returns   : a CoGe::Accessory::blastz_report::seq object
 Argument  : none
 Throws    :
 Comment   : the seq object has accessor functions for:
              file => file name
              start => start position
              stop => stop position
              length => length of sequence
              strand => strand used
              contig => something blastz uses when comparing contigs.  Not sure what this means, check blastz docs.

See Also   :

=cut

#################### subroutine header end ####################

sub seq2_info
  {
    my $self = shift;
    if (ref ($self->seq_info) =~ /hash/i)
      {
	return ($self->seq_info->{2});
      }
  }

sub parse_alignment
  {
    my $self = shift;
    my $data = shift;
    my $seg_count = 1;
    my $hsp_count =$self->hsp_count;
    $hsp_count++;
    my $hsp = new CoGe::Accessory::parse_report::HSP({
						      number=>$hsp_count,
						      query_name=>$self->seq_info->{1}->{header},
						      subject_name=>$self->seq_info->{2}->{header},
						      query_length=>$self->seq_info->{1}->{length},
						      subject_length=>$self->seq_info->{2}->{length},
						     });
    my @id;

    foreach my $line (split /\n/, $data)
      {
	next if $line =~ /{|}/; #skip brackets
	next unless $line;
	$line =~ s/^\s+//;
	my ($type, @data) = split /\s/, $line;
	$hsp->score($data[0]) if $type eq "s";
	if ($type eq "b")
	  {
	    $hsp->query_start($data[0]);
	    $hsp->subject_start($data[1]);
	  }
	if ($type eq "e")
	  {
	    $hsp->query_stop($data[0]);
	    $hsp->subject_stop($data[1]);
	  }
	if ($type eq "l")
	  {
	    my $qstart = $data[0];
	    my $qstop = $data[2];
	    my $sstart = $data[1];
	    my $sstop = $data[3];
 	    if ($self->seq1_info->strand =~ /-/)
 	      {
 		my $len = $self->seq1_info->length;
 		my $tmp = $qstart;
 		$qstart = ($len-$qstop+1);
 		$qstop  = ($len-$tmp+1);
 	      }
 	    elsif ($self->seq2_info->strand =~ /-/)
 	      {
 		my $len = $self->seq2_info->length;
 		my $tmp = $sstart;
 		$sstart = ($len-$sstop+1);
 		$sstop  = ($len-$tmp+1);
 	      }
	    my $seg = new CoGe::Accessory::blastz_report::segment({number=>$seg_count++,
								   query_start=>$qstart,
								   query_stop=>$qstop,
								   subject_start=>$sstart,
								   subject_stop=>$sstop,
								   identity=>$data[4],
								  });
	    push @id, [$data[2]-$data[0]+1, $data[4]];
	    push @{$hsp->segments},$seg;
	  }
      }
    my $strand = "++";
#    print STDERR Dumper $self->seq1_info, $self->seq2_info;
    if ($self->seq1_info->strand =~ /-/)
      {
	$strand = "+-";
	my $len = $self->seq1_info->length;
	my $tmp = $hsp->query_start;
	$hsp->query_start($len-$hsp->query_stop+1);
	$hsp->query_stop($len-$tmp+1);
      }
    elsif ($self->seq2_info->strand =~ /-/)
      {
	$strand = "+-";
	my $len = $self->seq2_info->length;
	my $tmp = $hsp->subject_start;
	$hsp->subject_start($len-$hsp->subject_stop+1);
	$hsp->subject_stop($len-$tmp+1);
      }
    $hsp->strand($strand);
    $hsp->eval("N/A");
    my $len = $hsp->query_stop-$hsp->query_start > $hsp->subject_stop-$hsp->subject_start ? $hsp->query_stop-$hsp->query_start+1 : $hsp->subject_stop-$hsp->subject_start+1;
    $hsp->length($len);
    my $pid =0;
    my $total=0;
    foreach my $item (@id)
      {
	$pid += ($item->[0]*$item->[1]);
	$total += $item->[0];
      }
    $hsp->match(ceil($pid/100));
    $hsp->percent_id(sprintf("%.2f",$pid/$total));
    $hsp->percent_sim($hsp->percent_id());
    $hsp->positive(ceil($pid/100));
    push @{$self->hsps},$hsp;
    $self->hsp_count($hsp_count);
  }

sub parse_sequences
  {
    my $self = shift;
    my $seq1 = $self->get_seq($self->seq1_info);
    my $seq2 = $self->get_seq($self->seq2_info);
    foreach my $hsp (@{$self->hsps})
      {
	my $align;
	my $seq_1;
	my $seq_2;
	my $qend;
	my $send;
	my $qgap=0;
	my $sgap=0;
	my $i =0;
	my $len = $#{$hsp->segments};
	foreach my $seg (sort {$a->qstart <=> $b->qstart} @{$hsp->segments})
	  {

	    my $length = $seg->query_stop-$seg->query_start+1;
	    $seg->length($length);
	    my $seg1 = substr($seq1, $seg->query_start-1, $length);
	    $seg->query_alignment($seg1);
	    my $seg2 = substr($seq2, $seg->subject_start-1, $length);
	    $seg->subject_alignment($seg2);
#	    print STDERR $seg->query_start,"-", $seg->query_end,"(",$seg->query_end- $seg->query_start,")\t";
#	    print STDERR $seg->subject_start,"-", $seg->subject_end,"(",$seg->subject_end- $seg->subject_start,")\n";
 	    if ($hsp->strand =~ /-/ && $seg2)
 	      {
# 		$seg2 = reverse($seg2);
# 		$seg2 =~ tr/ATCG/TAGC/;
 		if ($seq_2)
 		  {
 		    $seq_2 = $seg2.$seq_2;
 		  }
 		else
 		  {
 		    $seq_2 = $seg2;
 		  }
 	      }
 	    else
 	      {
		$seq_2 .= $seg2;
	      }
	    $send = $seg->subject_stop;
	    $seq_1 .= $seg1;
	    $qend = $seg->query_stop;
	    my $qi = $i;
	    my $qi1 = $i+1;
	    if ($hsp->strand =~ /-/)
	      {
		$qi = $len-$i;
		$qi1 = $qi-1;
	      }
#	    print STDERR "i: $i, qi: $qi, qi1: $qi1\n";
	    my $gap1 = $hsp->segments->[$i+1]->query_start-$hsp->segments->[$i]->query_end-1 if $hsp->segments->[$i+1];
	    my $gap2 = $hsp->segments->[$qi1]->subject_start-$hsp->segments->[$qi]->subject_end-1 if $hsp->segments->[$qi1] && $qi1 >= 0;
#	    print STDERR "\tqgap: ", $hsp->segments->[$i+1]->query_start,"-",$hsp->segments->[$i]->query_end,"::",$gap1,"\n" unless $i+1>$len;
#	    print STDERR "\tsgap: ", $hsp->segments->[$qi1]->subject_start,"-",$hsp->segments->[$qi]->subject_end,"::",$gap2,"\n" if $hsp->segments->[$qi1] && $qi1 >=0;

	    my ($qadd, $sadd);
	    if($gap1 && $i <= $len)
	      {
		$qadd = substr($seq1, $hsp->segments->[$i]->query_stop, $hsp->segments->[$i+1]->query_start-$hsp->segments->[$i]->query_end-1);
		$sadd = "-"x$gap1;
		$sgap += $gap1;
	      }
	    if ($gap2 && $qi1 >=0 )
	      {
		$sadd = substr($seq2, $hsp->segments->[$qi]->subject_stop, $hsp->segments->[$qi1]->subject_start-$hsp->segments->[$qi]->subject_end-1);
		$qadd = "-"x$gap2;
		$qgap += $gap2;
	      }

	    $seq_1 .= $qadd if $qadd;
 	    if ($hsp->strand =~ /-/ && $sadd)
 	      {

# 		$sadd = reverse($sadd);
# 		$sadd =~ tr/ATCG/TAGC/;
 		$seq_2 = $sadd.$seq_2;
 	      }
 	    else
 	      {
 		$seq_2 .= $sadd if $sadd;
 	      }
#	    print STDERR ($seq_1),"\n";
#	    print STDERR ( $seq_2),"\n";
	    $i++;
	  }
 	if ($hsp->strand =~ /-/)
 	  {
 	    $seq_2 = reverse($seq_2);
 	    $seq_2 =~ tr/ATCG/TAGC/;

 	  }
	my $length = length($seq_1) > length($seq_2) ? length($seq_1) : length($seq_2);
	foreach my $i (0..$length-1)
	  {
	    my $chr1 = substr($seq_1, $i, 1);
	    my $chr2 = substr ($seq_2, $i, 1);
	    exit unless ($chr1 && $chr2);
	    my $aln = ( $chr1 eq $chr2 ) ? "|" : " ";
	    $align .= $aln;
	  }
	    #$seg->alignment($aln);
	$hsp->alignment($align);
	$hsp->query_alignment($seq_1);
	$hsp->subject_alignment($seq_2);
	$hsp->query_gaps($qgap);
	$hsp->subject_gaps($sgap);
      }
  }

sub get_seq
  {
    my $self = shift;
    my $seq_obj = shift;
    return unless $seq_obj;
    my $file = $seq_obj->file;
    $file = $self->file_base.$file unless -r $file;
    unless (-r $file)
      {
	carp "can't read sequence file $file.  file_base=>".$self->file_base.", file=>".$seq_obj->file."\n";
      }
    my $seq;
    open (IN, $file);
    while (<IN>)
      {
	$seq .= $_;
      }
    close IN;
    my $name;
    ($name, $seq) = split/\n/, $seq ,2 if $seq =~ /^>/;
    $seq =~ s/\n//g;
    return $seq;
  }

#################### subroutine header begin ####################

=head2

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

1;

__END__
