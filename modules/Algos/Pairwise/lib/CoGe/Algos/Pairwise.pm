package CoGe::Algos::Pairwise;
use strict;
use base 'Class::Accessor';
use Data::Dumper;
use IO::Socket;
use CoGe::Accessory::Web;

BEGIN {
    use Exporter ();
    use vars qw($P $VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $NWALIGN $MATRIX_FILE);
    $VERSION     = '0.1';
    @ISA         = (@ISA, qw(Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
    __PACKAGE__->mk_accessors(qw(seqA seqB matrix gap gap_ext dpm alignA alignB nwalign nwalign_server_port));
#    $P = CoGe::Accessory::Web::get_defaults();
#    $NWALIGN = $P->{NWALIGN};
#    $MATRIX_FILE = $P->{BLASTMATRIX}."aa/BLOSUM62";
}

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module.
## You better edit it!

=head1 NAME

CoGe::Algos::Pairwise - Pairwise

=head1 SYNOPSIS

  use CoGe::Algos::Pairwise;
  my $seq1 = "PELICAN";
  my $seq2 = "CELLAICAN";
  my $pairwise = new CoGe::Algos::Pairwise
  $pairwise->seqA($seq1);
  $pairwise->seqB($seq2);

  #align the sequence
  my ($align1, $align2) = $pairwise->global_align();

  #pretty print the dynamic programming matrix
  $pairwise->print_dpm();

=head1 DESCRIPTION

Stub documentation for this module was created by ExtUtils::ModuleMaker.
It looks like the author of the extension was negligent enough
to leave the stub unedited.

Blah blah blah.

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

        Eric Lyons
	CPAN ID: MODAUTHOR
        UC Berkeley
	elyons(@t)nature.berkeley.edu

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

perl(1).

=cut

#################### main pod documentation end ###################

#################### subroutine header begin ####################

=head2 seqA

 Usage     : $pairwise->seqA($seq)
 Purpose   : get/set one of the sequences
 Returns   : string
 Argument  : string
 Throws    : none
 Comment   : this is how you set one of the two sequences to be aligned
           :

See Also   :

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 seqB

 Usage     : $pairwise->seqB($seq)
 Purpose   : get/set one of the sequences
 Returns   : string
 Argument  : string
 Throws    : none
 Comment   : this is how you set one of the two sequences to be aligned
           :

See Also   :

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 matrix

 Usage     : $pairwise->matrix($matrix)
 Purpose   : get/set the scoring matrix for the alignment
 Returns   : a ref to a hash of hash refs
 Argument  : a ref to a hash of hash refs
 Throws    : none
 Comment   : example matrix:  $matrix = {{A}=>{A=>1, T=>-1, C=>-1}, G=>-1},
           :                             {T}=>{A=>-1, T=>1, C=>-1}, G=>-1},
           :                             {C}=>{A=>-1, T=>-1, C=>1}, G=>-1},
           :                             {G}=>{A=>-1, T=>-1, C=>-1}, G=>1},};

           : if no matrix is supplied, this will use the BLOSUM62 by default

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 gap

 Usage     : $pairwise->gap(-10)
 Purpose   : get/set the gap opening cost
 Returns   : string/int
 Argument  : string/int
 Throws    : 0
 Comment   : If no gap is specified, -10 is used by default
           :

See Also   :

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 gap_ext

 Usage     : $pairwise->gap_ext(-2)
 Purpose   : get/set the gap extension cost
 Returns   : string/int
 Argument  : string/int
 Throws    : none
 Comment   : If no gap extension is specified, -2 is used by default
           :

See Also   :

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 dpm

 Usage     : my $dpm = $pairwise->dpm();
 Purpose   : storage place for the dynamic programming matrix used for the last alignment
 Returns   : a reference to a 2D matrix of hash refs
 Argument  : this is set by the global_align subroutine
 Throws    : none
 Comment   : If you want to get a hold the the DPM used to generate the alignment, this
           : is the puppy.  However, if you want to see it printed pretty, see the
           : print_dpm sub.
See Also   : sub print_dpm

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 alignA

 Usage     : $aligna = $pairwise->alignA();
 Purpose   : storage place for one of the aligned sequence post alignment
 Returns   : a string (or nothing if not set)
 Argument  : this is set by the alignment subroutine
 Throws    : none
 Comment   : Allows you to retrieve the alignment for one of the sequence
           : after the alignment has been run.

See Also   :

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 alignB

 Usage     : $alignb = $pairwise->alignB();
 Purpose   : storage place for one of the aligned sequence post alignment
 Returns   : a string (or nothing if not set)
 Argument  : this is set by the alignment subroutine
 Throws    : none
 Comment   : Allows you to retrieve the alignment for one of the sequence
           : after the alignment has been run.

See Also   :

=cut

#################### subroutine header end ####################

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

#################### subroutine header begin ####################

=head2 global_align

 Usage     : $pairwise->global_align()
 Purpose   : aligns two sequences stored in the pairwise object.  The sequences must have been previously
             set with $pairwise->seqA($seq) and $pairwise->seqB($seq2)
 Returns   : an array of two strings where each string is the global sequence alignment
 Argument  : None
 Throws    : returns 0 if either sequence is not defined
 Comment   : This does a global sequence alignment between two sequences using gap and gap extension
           : penalties.

See Also   :

=cut

#################### subroutine header end ####################

sub global_align {
    my $self = shift;
    my %opts = @_;

    my $seq1 = $opts{seqA};
    $seq1 = $self->seqA unless defined $seq1;
    my $seq2 = $opts{seqB};
    $seq2 = $self->seqB unless defined $seq2;
    return 0 unless ( $seq1 && $seq2 );
    my $gap = $opts{gap};
    $gap = $self->gap unless defined $gap;
    my $gap_ext = $opts{gap_ext};
    $gap_ext = $self->gap_ext unless defined $gap_ext;

    $gap = -10 unless defined $gap;
    $gap_ext = -2 unless $gap_ext;
    my $matrix = $opts{matrix};  #path to blast formated alignment matrix;

    $P = CoGe::Accessory::Web::get_defaults($opts{config});
    $NWALIGN = $P->{NWALIGN};
    $MATRIX_FILE = $P->{BLASTMATRIX}."aa/BLOSUM62";
    $matrix = $MATRIX_FILE unless $matrix && -r $matrix;

    my ($align1, $align2);
    if ($self->nwalign_server_port)
      {
	my $sock = IO::Socket::INET->new(
					 PeerAddr => 'localhost:'.$self->nwalign_server_port,
					) or die "Can't bind: $@\n";
	my $cmd = " --matrix $matrix --gap_extend $gap_ext --gap_open $gap $seq1 $seq2";
	print $sock $cmd;
	my $res;
	$sock->recv($res, 128000);
	($align1, $align2) = split/\s+/, $res, 2;
      }
    else
      {
	my $prog = $self->nwalign;
	$prog = $NWALIGN unless $prog;
	my $cmd = "$prog --matrix=$matrix --gap_extend=$gap_ext --gap_open=$gap $seq1 $seq2";
	open (RUN, "$cmd |");
	($align1, $align2) = <RUN>;
	close RUN;
      }
    $align1 =~ s/\n//;
    $align2 =~ s/\n//;
    return ( $align1, $align2 );
}

sub global_align_perl {
    ########################
    # initialization stage #
    ########################
    ## set up an empty matrix, fill is with our starting values
    ## make horizontals point "left", and verticals point up
    ## row or column value is the gap penalty multiplied by the position in
    ## the matrix
    my ($self) = (@_);
    my $seq1   = $self->seqA;
    my $seq2   = $self->seqB;
    return 0 unless ( $seq1 && $seq2 );
    $self->gap(-10)    unless $self->gap;
    $self->gap_ext(-2) unless $self->gap_ext;
    my $GAP     = $self->gap;
    my $GAP_EXT = $self->gap_ext;
    $self->_initialize_default_scoring_matrix unless $self->matrix;
    my $sm = $self->matrix;

    my @matrix;
    $matrix[0][0]{score}   = 0;
    $matrix[0][0]{pointer} = "none";
    $matrix[0][0]{gap}     = 0;
    for ( my $j = 1 ; $j <= length($seq1) ; $j++ ) {
        my $cost = $matrix[0][ $j - 1 ]{gap} ? $GAP_EXT : $GAP;
        $matrix[0][$j]{score}   = $matrix[0][ $j - 1 ]{score} + $cost;
        $matrix[0][$j]{pointer} = "lt";
        $matrix[0][$j]{gap}     = 1;
    }
    for ( my $i = 1 ; $i <= length($seq2) ; $i++ ) {
        my $cost = $matrix[ $i - 1 ][0]{gap} ? $GAP_EXT : $GAP;
        $matrix[$i][0]{score}   = $matrix[ $i - 1 ][0]{score} + $cost;
        $matrix[$i][0]{pointer} = "up";
        $matrix[$i][0]{gap}     = 1;
    }

    ########################
    # fill stage           #
    ########################
    ## run through every element and calculate whether or not there is a
    ## match, score accordingly
    ##

    for ( my $i = 1 ; $i <= length($seq2) ; $i++ ) {
        for ( my $j = 1 ; $j <= length($seq1) ; $j++ ) {
            my ( $diagonal_score, $left_score, $up_score );

            # calculate match score
            my $letter1 = substr( $seq1, $j - 1, 1 );
            my $letter2 = substr( $seq2, $i - 1, 1 );
            my $match_score = $sm->{$letter1}{$letter2};
            $match_score = 0
              unless defined $match_score
            ;    #some odd character?  Let's not make it hurt or add
            $diagonal_score = $matrix[ $i - 1 ][ $j - 1 ]{score} + $match_score;

            # calculate gap scores
            my $cost = $matrix[ $i - 1 ][$j]{gap} ? $GAP_EXT : $GAP;
            $up_score   = $matrix[ $i - 1 ][$j]{score} + $cost;
            $cost       = $matrix[$i][ $j - 1 ]{gap} ? $GAP_EXT : $GAP;
            $left_score = $matrix[$i][ $j - 1 ]{score} + $cost;

            # choose best score
            if ( $diagonal_score >= $up_score ) {
                if ( $diagonal_score >= $left_score ) {
                    $matrix[$i][$j]{score}   = $diagonal_score;
                    $matrix[$i][$j]{pointer} = "dg";
                    $matrix[$i][$j]{gap}     = 0;
                }
                else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "lt";
                    $matrix[$i][$j]{gap}     = 1;
                }
            }
            else {
                if ( $up_score >= $left_score ) {
                    $matrix[$i][$j]{score}   = $up_score;
                    $matrix[$i][$j]{pointer} = "up";
                    $matrix[$i][$j]{gap}     = 1;
                }
                else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "lt";
                    $matrix[$i][$j]{gap}     = 1;
                }
            }
        }
    }

    ######################
    # trace-back         #
    ######################

    my $align1 = "";    # reserve some global variables here
    my $align2 = "";

    # start at last cell of matrix, (i,j)
    my $j = length($seq1);
    my $i = length($seq2);

    ## while loop condition is always true
    while (1) {
      last if $matrix[$i][$j]{pointer} eq "none"; # exits while loop
      # at first cell of matrix

      if ($matrix[$i][$j]{pointer} eq "dg") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= substr($seq2, $i-1, 1);
        $i--;			#decrement operator
        $j--;			#decrement operator
      } elsif ($matrix[$i][$j]{pointer} eq "lt") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= "-";
        $j--;			#decrement operator
      } elsif ($matrix[$i][$j]{pointer} eq "up") {
        $align1 .= "-";
        $align2 .= substr($seq2, $i-1, 1);
        $i--;			#decrement operator
      }
    }
    $self->dpm( \@matrix );
    $align1 = reverse($align1);
    $align2 = reverse($align2);
    $self->alignA($align1);
    $self->alignB($align2);
    return ( $align1, $align2 );
}

#################### subroutine header begin ####################

=head2 print_dpm

 Usage     : $self->print_dpm
 Purpose   : pretty prints the dynamic programming matrix
 Returns   : none
 Argument  : none
 Throws    : 0
 Comment   : output is:  score:gap_flag:trace_direction
           : e.g.  -10:1:up         -3:0:dg        -13:1:lt
           : lt = trace is from left cell
           : up = trace is from above cell
           : gd = trace is from above diagonal cell

See Also   :

=cut

#################### subroutine header end ####################

sub print_dpm {
    my $self = shift;
    my $matrix = shift || $self->dpm();
    return 0 unless $matrix;
    foreach my $i (@$matrix) {
        foreach my $j (@$i) {
            print join( ":",
                sprintf( "%4d", $j->{score} ),
                $j->{gap}, $j->{pointer} );
            print "\t";
        }
        print "\n";
    }
}

#################### subroutine header begin ####################

=head2 print_align

 Usage     : $pw->print_align();
 Purpose   : prints a pretty alignment
 Returns   : none
 Argument  : string (int) for the number of characters before wrapping the
           : alignment to the next string
 Throws    : none
 Comment   : a simple way to get a pretty and easy to read alignment
           :

See Also   :

=cut

#################### subroutine header end ####################

sub print_align {
    my $self = shift;
    my %opts = @_;
    my $wrap = $opts{wrap} || 60;
    my $seqA = $opts{seqA} || $opts{seq1} || $self->alignA;
    my $seqB = $opts{seqB} || $opts{seq2} || $self->alignB;
    my @seqA = split //, $seqA;
    my @seqB = split //, $seqB;
    my ( $s1, $s2, $s3 );
    my $count = 1;
    my $index = 0;

    foreach my $c1 (@seqA) {
        my $c2 = $seqB[$index];
        $s1 .= $c1;
        $s2 .= $c2;
        if ( $c1 eq $c2 ) {
            $s3 .= ":";
        }
        elsif ($self->matrix->{$c1}
            && defined $self->matrix->{$c1}{$c2}
            && $self->matrix->{$c1}{$c2} >= 0 )
        {
            $s3 .= ".";
        }
        else {
            $s3 .= " ";
        }
        $count++;
        $index++;
        if ( $count == $wrap ) {
            print join( "\n", $s1, $s2, $s3, "\n" );
            $count = 1;
            $s1    = undef;
            $s2    = undef;
            $s3    = undef;
        }
    }
    print join( "\n", $s1, $s2, $s3, "\n" );
}

#################### subroutine header begin ####################

=head2 _initialize_default_scoring_matrix

 Usage     : $pairwise->_initialize_default_scoring_matrix
 Purpose   : set the scoring matrix to BLOSOM62
 Returns   : $self->matrix()
 Argument  : none
 Throws    : none
 Comment   : if no scoring matrix has been specified, this matrix is used
           :

See Also   :

=cut

#################### subroutine header end ####################

sub _initialize_default_scoring_matrix {
    #BLOSUM62
    my $self = shift;
    $self->matrix(
        {
            'S' => {
                'S' => '4',
                'F' => '-2',
                'T' => '1',
                'N' => '1',
                'K' => '0',
                '*' => '-4',
                'Y' => '-2',
                'E' => '0',
                'V' => '-2',
                'Z' => '0',
                'Q' => '0',
                'M' => '-1',
                'C' => '-1',
                'L' => '-2',
                'A' => '1',
                'W' => '-3',
                'X' => '0',
                'P' => '-1',
                'B' => '0',
                'H' => '-1',
                'D' => '0',
                'R' => '-1',
                'I' => '-2',
                'G' => '0'
            },
            'F' => {
                'S' => '-2',
                'F' => '6',
                'T' => '-2',
                'N' => '-3',
                'K' => '-3',
                '*' => '-4',
                'Y' => '3',
                'E' => '-3',
                'V' => '-1',
                'Z' => '-3',
                'Q' => '-3',
                'M' => '0',
                'C' => '-2',
                'L' => '0',
                'A' => '-2',
                'W' => '1',
                'X' => '-1',
                'P' => '-4',
                'B' => '-3',
                'H' => '-1',
                'D' => '-3',
                'R' => '-3',
                'I' => '0',
                'G' => '-3'
            },
            'T' => {
                'S' => '1',
                'F' => '-2',
                'T' => '5',
                'N' => '0',
                'K' => '-1',
                '*' => '-4',
                'Y' => '-2',
                'E' => '-1',
                'V' => '0',
                'Z' => '-1',
                'Q' => '-1',
                'M' => '-1',
                'C' => '-1',
                'L' => '-1',
                'A' => '0',
                'W' => '-2',
                'X' => '0',
                'P' => '-1',
                'B' => '-1',
                'H' => '-2',
                'D' => '-1',
                'R' => '-1',
                'I' => '-1',
                'G' => '-2'
            },
            'N' => {
                'S' => '1',
                'F' => '-3',
                'T' => '0',
                'N' => '6',
                'K' => '0',
                '*' => '-4',
                'Y' => '-2',
                'E' => '0',
                'V' => '-3',
                'Z' => '0',
                'Q' => '0',
                'M' => '-2',
                'C' => '-3',
                'L' => '-3',
                'A' => '-2',
                'W' => '-4',
                'X' => '-1',
                'P' => '-2',
                'B' => '3',
                'H' => '1',
                'D' => '1',
                'R' => '0',
                'I' => '-3',
                'G' => '0'
            },
            'K' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-1',
                'N' => '0',
                'K' => '5',
                '*' => '-4',
                'Y' => '-2',
                'E' => '1',
                'V' => '-2',
                'Z' => '1',
                'Q' => '1',
                'M' => '-1',
                'C' => '-3',
                'L' => '-2',
                'A' => '-1',
                'W' => '-3',
                'X' => '-1',
                'P' => '-1',
                'B' => '0',
                'H' => '-1',
                'D' => '-1',
                'R' => '2',
                'I' => '-3',
                'G' => '-2'
            },
            '*' => {
                'S' => '-4',
                'F' => '-4',
                'T' => '-4',
                'N' => '-4',
                'K' => '-4',
                '*' => '1',
                'Y' => '-4',
                'E' => '-4',
                'V' => '-4',
                'Z' => '-4',
                'Q' => '-4',
                'M' => '-4',
                'C' => '-4',
                'L' => '-4',
                'A' => '-4',
                'W' => '-4',
                'X' => '-4',
                'P' => '-4',
                'B' => '-4',
                'H' => '-4',
                'D' => '-4',
                'R' => '-4',
                'I' => '-4',
                'G' => '-4'
            },
            'Y' => {
                'S' => '-2',
                'F' => '3',
                'T' => '-2',
                'N' => '-2',
                'K' => '-2',
                '*' => '-4',
                'Y' => '7',
                'E' => '-2',
                'V' => '-1',
                'Z' => '-2',
                'Q' => '-1',
                'M' => '-1',
                'C' => '-2',
                'L' => '-1',
                'A' => '-2',
                'W' => '2',
                'X' => '-1',
                'P' => '-3',
                'B' => '-3',
                'H' => '2',
                'D' => '-3',
                'R' => '-2',
                'I' => '-1',
                'G' => '-3'
            },
            'E' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-1',
                'N' => '0',
                'K' => '1',
                '*' => '-4',
                'Y' => '-2',
                'E' => '5',
                'V' => '-2',
                'Z' => '4',
                'Q' => '2',
                'M' => '-2',
                'C' => '-4',
                'L' => '-3',
                'A' => '-1',
                'W' => '-3',
                'X' => '-1',
                'P' => '-1',
                'B' => '1',
                'H' => '0',
                'D' => '2',
                'R' => '0',
                'I' => '-3',
                'G' => '-2'
            },
            'V' => {
                'S' => '-2',
                'F' => '-1',
                'T' => '0',
                'N' => '-3',
                'K' => '-2',
                '*' => '-4',
                'Y' => '-1',
                'E' => '-2',
                'V' => '4',
                'Z' => '-2',
                'Q' => '-2',
                'M' => '1',
                'C' => '-1',
                'L' => '1',
                'A' => '0',
                'W' => '-3',
                'X' => '-1',
                'P' => '-2',
                'B' => '-3',
                'H' => '-3',
                'D' => '-3',
                'R' => '-3',
                'I' => '3',
                'G' => '-3'
            },
            'Z' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-1',
                'N' => '0',
                'K' => '1',
                '*' => '-4',
                'Y' => '-2',
                'E' => '4',
                'V' => '-2',
                'Z' => '4',
                'Q' => '3',
                'M' => '-1',
                'C' => '-3',
                'L' => '-3',
                'A' => '-1',
                'W' => '-3',
                'X' => '-1',
                'P' => '-1',
                'B' => '1',
                'H' => '0',
                'D' => '1',
                'R' => '0',
                'I' => '-3',
                'G' => '-2'
            },
            'Q' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-1',
                'N' => '0',
                'K' => '1',
                '*' => '-4',
                'Y' => '-1',
                'E' => '2',
                'V' => '-2',
                'Z' => '3',
                'Q' => '5',
                'M' => '0',
                'C' => '-3',
                'L' => '-2',
                'A' => '-1',
                'W' => '-2',
                'X' => '-1',
                'P' => '-1',
                'B' => '0',
                'H' => '0',
                'D' => '0',
                'R' => '1',
                'I' => '-3',
                'G' => '-2'
            },
            'M' => {
                'S' => '-1',
                'F' => '0',
                'T' => '-1',
                'N' => '-2',
                'K' => '-1',
                '*' => '-4',
                'Y' => '-1',
                'E' => '-2',
                'V' => '1',
                'Z' => '-1',
                'Q' => '0',
                'M' => '5',
                'C' => '-1',
                'L' => '2',
                'A' => '-1',
                'W' => '-1',
                'X' => '-1',
                'P' => '-2',
                'B' => '-3',
                'H' => '-2',
                'D' => '-3',
                'R' => '-1',
                'I' => '1',
                'G' => '-3'
            },
            'C' => {
                'S' => '-1',
                'F' => '-2',
                'T' => '-1',
                'N' => '-3',
                'K' => '-3',
                '*' => '-4',
                'Y' => '-2',
                'E' => '-4',
                'V' => '-1',
                'Z' => '-3',
                'Q' => '-3',
                'M' => '-1',
                'C' => '9',
                'L' => '-1',
                'A' => '0',
                'W' => '-2',
                'X' => '-2',
                'P' => '-3',
                'B' => '-3',
                'H' => '-3',
                'D' => '-3',
                'R' => '-3',
                'I' => '-1',
                'G' => '-3'
            },
            'L' => {
                'S' => '-2',
                'F' => '0',
                'T' => '-1',
                'N' => '-3',
                'K' => '-2',
                '*' => '-4',
                'Y' => '-1',
                'E' => '-3',
                'V' => '1',
                'Z' => '-3',
                'Q' => '-2',
                'M' => '2',
                'C' => '-1',
                'L' => '4',
                'A' => '-1',
                'W' => '-2',
                'X' => '-1',
                'P' => '-3',
                'B' => '-4',
                'H' => '-3',
                'D' => '-4',
                'R' => '-2',
                'I' => '2',
                'G' => '-4'
            },
            'A' => {
                'S' => '1',
                'F' => '-2',
                'T' => '0',
                'N' => '-2',
                'K' => '-1',
                '*' => '-4',
                'Y' => '-2',
                'E' => '-1',
                'V' => '0',
                'Z' => '-1',
                'Q' => '-1',
                'M' => '-1',
                'C' => '0',
                'L' => '-1',
                'A' => '4',
                'W' => '-3',
                'X' => '0',
                'P' => '-1',
                'B' => '-2',
                'H' => '-2',
                'D' => '-2',
                'R' => '-1',
                'I' => '-1',
                'G' => '0'
            },
            'W' => {
                'S' => '-3',
                'F' => '1',
                'T' => '-2',
                'N' => '-4',
                'K' => '-3',
                '*' => '-4',
                'Y' => '2',
                'E' => '-3',
                'V' => '-3',
                'Z' => '-3',
                'Q' => '-2',
                'M' => '-1',
                'C' => '-2',
                'L' => '-2',
                'A' => '-3',
                'W' => '11',
                'X' => '-2',
                'P' => '-4',
                'B' => '-4',
                'H' => '-2',
                'D' => '-4',
                'R' => '-3',
                'I' => '-3',
                'G' => '-2'
            },
            'X' => {
                'S' => '0',
                'F' => '-1',
                'T' => '0',
                'N' => '-1',
                'K' => '-1',
                '*' => '-4',
                'Y' => '-1',
                'E' => '-1',
                'V' => '-1',
                'Z' => '-1',
                'Q' => '-1',
                'M' => '-1',
                'C' => '-2',
                'L' => '-1',
                'A' => '0',
                'W' => '-2',
                'X' => '-1',
                'P' => '-2',
                'B' => '-1',
                'H' => '-1',
                'D' => '-1',
                'R' => '-1',
                'I' => '-1',
                'G' => '-1'
            },
            'P' => {
                'S' => '-1',
                'F' => '-4',
                'T' => '-1',
                'N' => '-2',
                'K' => '-1',
                '*' => '-4',
                'Y' => '-3',
                'E' => '-1',
                'V' => '-2',
                'Z' => '-1',
                'Q' => '-1',
                'M' => '-2',
                'C' => '-3',
                'L' => '-3',
                'A' => '-1',
                'W' => '-4',
                'X' => '-2',
                'P' => '7',
                'B' => '-2',
                'H' => '-2',
                'D' => '-1',
                'R' => '-2',
                'I' => '-3',
                'G' => '-2'
            },
            'B' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-1',
                'N' => '3',
                'K' => '0',
                '*' => '-4',
                'Y' => '-3',
                'E' => '1',
                'V' => '-3',
                'Z' => '1',
                'Q' => '0',
                'M' => '-3',
                'C' => '-3',
                'L' => '-4',
                'A' => '-2',
                'W' => '-4',
                'X' => '-1',
                'P' => '-2',
                'B' => '4',
                'H' => '0',
                'D' => '4',
                'R' => '-1',
                'I' => '-3',
                'G' => '-1'
            },
            'H' => {
                'S' => '-1',
                'F' => '-1',
                'T' => '-2',
                'N' => '1',
                'K' => '-1',
                '*' => '-4',
                'Y' => '2',
                'E' => '0',
                'V' => '-3',
                'Z' => '0',
                'Q' => '0',
                'M' => '-2',
                'C' => '-3',
                'L' => '-3',
                'A' => '-2',
                'W' => '-2',
                'X' => '-1',
                'P' => '-2',
                'B' => '0',
                'H' => '8',
                'D' => '-1',
                'R' => '0',
                'I' => '-3',
                'G' => '-2'
            },
            'D' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-1',
                'N' => '1',
                'K' => '-1',
                '*' => '-4',
                'Y' => '-3',
                'E' => '2',
                'V' => '-3',
                'Z' => '1',
                'Q' => '0',
                'M' => '-3',
                'C' => '-3',
                'L' => '-4',
                'A' => '-2',
                'W' => '-4',
                'X' => '-1',
                'P' => '-1',
                'B' => '4',
                'H' => '-1',
                'D' => '6',
                'R' => '-2',
                'I' => '-3',
                'G' => '-1'
            },
            'R' => {
                'S' => '-1',
                'F' => '-3',
                'T' => '-1',
                'N' => '0',
                'K' => '2',
                '*' => '-4',
                'Y' => '-2',
                'E' => '0',
                'V' => '-3',
                'Z' => '0',
                'Q' => '1',
                'M' => '-1',
                'C' => '-3',
                'L' => '-2',
                'A' => '-1',
                'W' => '-3',
                'X' => '-1',
                'P' => '-2',
                'B' => '-1',
                'H' => '0',
                'D' => '-2',
                'R' => '5',
                'I' => '-3',
                'G' => '-2'
            },
            'I' => {
                'S' => '-2',
                'F' => '0',
                'T' => '-1',
                'N' => '-3',
                'K' => '-3',
                '*' => '-4',
                'Y' => '-1',
                'E' => '-3',
                'V' => '3',
                'Z' => '-3',
                'Q' => '-3',
                'M' => '1',
                'C' => '-1',
                'L' => '2',
                'A' => '-1',
                'W' => '-3',
                'X' => '-1',
                'P' => '-3',
                'B' => '-3',
                'H' => '-3',
                'D' => '-3',
                'R' => '-3',
                'I' => '4',
                'G' => '-4'
            },
            'G' => {
                'S' => '0',
                'F' => '-3',
                'T' => '-2',
                'N' => '0',
                'K' => '-2',
                '*' => '-4',
                'Y' => '-3',
                'E' => '-2',
                'V' => '-3',
                'Z' => '-2',
                'Q' => '-2',
                'M' => '-3',
                'C' => '-3',
                'L' => '-4',
                'A' => '0',
                'W' => '-2',
                'X' => '-1',
                'P' => '-2',
                'B' => '-1',
                'H' => '-2',
                'D' => '-1',
                'R' => '-2',
                'I' => '-4',
                'G' => '6'
            }
        }
    );
    return $self->matrix;
}

1;
# The preceding line will help the module return a true value
