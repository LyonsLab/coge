package CoGe::Algos::KsCalc;
use strict;
use base qw(CoGe::Algos::Pairwise);
#use CoGeX;
use CoGe::Algos::Codeml;
use Data::Dumper;
use Carp;
use Benchmark qw(:all);

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION = '0.1';
    @ISA = ( @ISA, qw(Exporter) );
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
    __PACKAGE__->mk_accessors(qw(version name1 name2 prot1 prot2 palign1 palign2 dna1 dna2 dalign1 dalign2 gaplessP1 gaplessP2 gaplessD1 gaplessD2 results gapless_prot_pid gapless_DNA_pid prot_pid DNA_pid feat1 feat2 benchmark tmpdir debug config));
}

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module.
## You better edit it!

=head1 NAME

CoGe::Algos::KsCalc - CoGe::Algos::KsCalc

=head1 SYNOPSIS

  use CoGe::Algos::KsCalc;

  my $ks = CoGe::Algos::KsCalc->new ();

  $ks->feat1($coge_feat1);  #a CoGeX::Feature object of type CDS
  $ks->feat2($coge_feat2);  #a CoGeX::Feature object of type CDS

  my $res = $ks->KsCalc("seq.align");
  print "Ka = ", $res->{'dN'},"\n";
  print "Ks = ", $res->{'dS'},"\n";
  print "Ka/Ks = ", $res->{'dN/dS'},"\n";

=head1 DESCRIPTION

Inherets from CoGe::Algos::Pairwise and provides extended funcitonality
to calculate Ks Kd (synonymous and nonsynonymous substitution rates) from sequences
in the CoGeX database.  Calculations are performed by Codeml of PAML
(Phylogenetic Analysis by Maximum Likelihood) package of Ziheng Yang.  See
http://abacus.gene.ucl.ac.uk/software/paml.html for more information.

This object is to make it easy to calculate Ks and Ka values from the genomes database
for any feature that has a protein sequence and a name.

It procedure is:
1. for each sequence name, find the longest protein sequence assoicate with that name
2. align those sequences using a global pairwise alignment algorithm (default parameters
   set in CoGe::Algos::Pairwise)
3. remove gaps from alignment
4. generate corresponding DNA sequence from gapless protein alignment
5. run Codeml using Comp_Geomics::Algos::Codeml
6. store results in $self->results as a hash returned from CoGe::Algos::Codeml
   (see example above)

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

=head2 new

 Usage     : my $ks = new CoGe::Algos::KsCalc
 Purpose   : creates a KsCalc object
 Returns   : KsCalc object
 Argument  : none
 Throws    : none
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub new {
    my ( $class, %parameters ) = @_;

    my $self = bless( {}, ref($class) || $class );
    $self->init();
    return $self;
}

=head2 Class::Accessor functions

This is a list of the Class::Accessor functions and the information they hold:
version          Version of the data source to limit search when finding sequences
                 for sequence names.  This is in the data_information table of
                 the genomes database
name1            The name of one of the sequences in the genomes database for which
                 a sequence is to be searched
name2            The name of one of the sequences in the genomes database for which
                 a sequence is to be searched
prot1            Storage for the protein sequence for name1
prot2            Storage for the protein sequence for name2
palign1          Storage for the protein alignment for prot1 after alignment
palign2          Storage for the protein alignment for prot2 after alignment
dna1             Storage for the dna sequence for name1
dna2             Storage for the dna sequence for name2
gaplessP1        Storage for the gapless alignment of palign1
gaplessP2        Storage for the gapless alignment of palign2
gaplessD1        Storage for the gapless alignment of dna1 based on gaplessP1
gaplessD2        Storage for the gapless alignment of dna2 based on gaplessP2
results          Storage for the results hash for the Ks calculation by codeml

feat1            Storage for CoGeX feature object for feature 1
feat2            Storage for CoGeX feature object for feature 2
=cut

#################### subroutine header begin ####################

=head2 init

 Usage     : $self->init (called by new)
 Purpose   : create and set default parameters.  Currently, sets tmpdir to /tmp/

 Returns   : none
 Argument  : none
 Throws    : none
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub init {
    my $self = shift;
    $self->tmpdir("/tmp/");
}

#################### subroutine header begin ####################

=head2 palign

 Usage     : $ks->palign
 Purpose   : after the names have been set, generates the global alignment
             ($self->global_align inhereted from Pairwise) and saves the aligned
             sequences in $self->palign1 and $self->palign2
 Returns   : 1 if successful, 0 if not
 Argument  : none
 Throws    : 0
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub palign {
    my $self = shift;
    return 0 unless $self->_check_seqs();

    $self->seqA($self->prot1);
    $self->seqB($self->prot2);
#    print STDERR "KsCalc:  ".$self->config,"\n";
    my ($align1, $align2) = $self->global_align(config=>$self->config);
    $self->palign1($align1);
    $self->palign2($align2);
    $self->_generate_DNA_alignment();
    $self->_generate_gapless();
    $self->gapless_DNA_calc_pid();
    $self->gapless_prot_calc_pid();
    $self->DNA_calc_pid();
    $self->prot_calc_pid();
    return 1 if ( $self->palign1() && $self->palign2() );
    return 0;
}

#################### subroutine header begin ####################

=head2 gapless_DNA_calc_pid

 Usage     : my $pid = $ks->gapless_DNA_calc_pid();
 Purpose   : Calculates the percent identity between the two stored aligned gapless nucleotide aligned sequences
 Returns   : pid a number between 0 and 100, no cropping of significant figures
 Argument  : none
 Throws    : undef + errors to STDOUT if either sequence is not defined
 Comment   : calls $self->calc_pid
           :
See Also   : calc_pid

=cut

#################### subroutine header end ####################

sub gapless_DNA_calc_pid {
    my $self = shift;
    my ( $seq1, $seq2 ) = ( $self->gaplessD1, $self->gaplessD2 );
    unless ( $seq1 && $seq2 ) {
        print STDERR
'Missing at least one of the gapless DNA sequences.  Perhaps $self->palign was not run?';
        return;
    }
    $self->gapless_DNA_pid( $self->calc_pid( $seq1, $seq2 ) );
    return ( $self->gapless_DNA_pid );
}

#################### subroutine header begin ####################

=head2 gapless_prot_calc_pid

 Usage     : my $pid = $ks->gapless_prot_calc_pid();
 Purpose   : Calculates the percent identity between the two stored aligned gapless protein aligned sequences
 Returns   : pid a number between 0 and 100, no cropping of significant figures
 Argument  : none
 Throws    : undef + errors to STDOUT if either sequence is not defined
 Comment   : calls $self->calc_pid
           :
See Also   : calc_pid

=cut

#################### subroutine header end ####################

sub gapless_prot_calc_pid {
    my $self = shift;
    my ( $seq1, $seq2 ) = ( $self->gaplessP1, $self->gaplessP2 );
    unless ( $seq1 && $seq2 ) {
        print STDERR
'Missing at least one of the gapless protein sequences.  Perhaps $self->palign was not run?';
        return;
    }
    $self->gapless_prot_pid( $self->calc_pid( $seq1, $seq2 ) );
    return ( $self->gapless_prot_pid );
}

#################### subroutine header begin ####################

=head2 DNA_calc_pid -- NOT IMPLEMENTED!

 Usage     : my $pid = $ks->DNA_calc_pid();
 Purpose   : Calculates the percent identity between the two stored aligned nucleotide aligned sequences
 Returns   : pid a number between 0 and 100, no cropping of significant figures
 Argument  : none
 Throws    : undef + errors to STDOUT if either sequence is not defined
 Comment   : calls $self->calc_pid
           :
See Also   : calc_pid

=cut

#################### subroutine header end ####################

sub DNA_calc_pid {
    my $self = shift;
    my ( $seq1, $seq2 ) = ( $self->dalign1, $self->dalign2 );
    unless ( $seq1 && $seq2 ) {
        print STDERR
'Missing at least one of the DNA sequences.  Perhaps $self->palign was not run?';
        return;
    }
    $self->DNA_pid( $self->calc_pid( $seq1, $seq2 ) );
    return ( $self->DNA_pid );
}

#################### subroutine header begin ####################

=head2 prot_calc_pid

 Usage     : my $pid = $ks->prot_calc_pid();
 Purpose   : Calculates the percent identity between the two stored aligned protein aligned sequences
 Returns   : pid a number between 0 and 100, no cropping of significant figures
 Argument  : none
 Throws    : undef + errors to STDOUT if either sequence is not defined
 Comment   : calls $self->calc_pid
           :
See Also   : calc_pid

=cut

#################### subroutine header end ####################

sub prot_calc_pid {
    my $self = shift;
    my ( $seq1, $seq2 ) = ( $self->palign1, $self->palign2 );
    unless ( $seq1 && $seq2 ) {
        print STDERR
'Missing at least one of the protein sequences.  Perhaps $self->palign was not run?';
        return;
    }
    $self->prot_pid( $self->calc_pid( $seq1, $seq2 ) );
    return ( $self->prot_pid );
}

#################### subroutine header begin ####################

=head2 calc_pid

 Usage     : my $pid = $ks->calc_pid($seq1, $seq2);
 Purpose   : Calculates the percent identity between two sequence
 Returns   : pid a number between 0 and 100, no cropping of significant figures
 Argument  : two strings (two sequences)
 Throws    : undef + errors to STDOUT if either sequence is not defined or
             if the sequences are not of equal length
 Comment   : The calculation of Num_identical_charaters / total_Num_characters
           : uses all non-gap characters from the $seq1 to calculate total_Num_characters
See Also   :

=cut

#################### subroutine header end ####################

sub calc_pid {
    my $self = shift;
    my ( $seq1, $seq2 ) = @_;
    unless ( $seq1 && $seq2 ) {
        carp "Sequences not valid in calc_pid";
        return;
    }
    unless ( length($seq1) == length($seq2) ) {
        print STDERR
"sequences are of different lengths.  Percent identity calculation may be incorrect.\n"
          if $self->debug;
    #	print STDERR join ("\t", $self->feat1->names),"\n" if $self->feat1->names;
    #	print STDERR $seq1,"\n\n";
    #	print STDERR join ("\t", $self->feat2->names),"\n" if $self->feat2->names;
    #	print STDERR $seq2,"\n\n"
    }
    my @seq1 = split //, $seq1;
    my @seq2 = split //, $seq2;
    my $total = 0;
    my $id    = 0;
    for ( my $i = 0 ; $i < scalar @seq1 ; $i++ ) {
        $total++ unless $seq1[$i] eq "-";
        next unless $seq2[$i];
        $id++ if $seq1[$i] eq $seq2[$i];
    }
    return $id / $total * 100;
}

#################### subroutine header begin ####################

=head2 KsCalc

 Usage     : $self->KsCalc
 Purpose   : runs $self->palign and runs Codeml from
             CoGe::Algos::Codeml and saves the results in
             $self->results
 Returns   : hash ref of results:
             'dN/dS' => non-synonymous over sysnonymous substitution
             'dN'    => non-synonymous substitution
             'dS'    => synonymous substitution
             'pID'   => percent identical
 Argument  : none
 Throws    : 0 if there was a problem running alignment
 Comment   : This is the mama-jama of this module
           :

See Also   :

=cut

#################### subroutine header end ####################

sub KsCalc
  {
    my $self = shift;
    my %opts = @_;
    my $config = $opts{config};
    $config = $self->config unless $config;
    $self->config($config) if $config; #set the self environmental var if set;
    my $t1 = new Benchmark if $self->benchmark;
    my $tmp = $self->palign(); #returns 0 if fails
    my $t2 = new Benchmark if $self->benchmark;
    unless ($tmp)
      {
	carp ("Problem running alignment");
	return 0;
      }
    my $cml = new CoGe::Algos::Codeml(alignment=>$self->phylip_align, config=>$config);
    $cml->run();
    my $t3 = new Benchmark if $self->benchmark;
    $self->results( $cml->results );
    if ( $self->benchmark ) {
        my $align_time  = timestr( timediff( $t2, $t1 ) );
        my $codeml_time = timestr( timediff( $t3, $t2 ) );
        print qq{
Time to align:     $align_time
Time to codeml:    $codeml_time
};
    }

    return $self->results;

}

#################### subroutine header begin ####################

=head2 phylip_align

 Usage     : $self->phylip_align('dna')
 Purpose   : generates a phylip format alignment of the gapless sequences
             default is to generate the dna gapless alignment, but
             you can specify 'prot' as an argument and the gapless protein
             alignment will be generated
 Returns   : a string
 Argument  : optional ('prot' | 'dna') default is 'dna'
 Throws    : none
 Comment   : this is used to generate the alignment file that Codeml uses
           : for its Ks calculation

See Also   :

=cut

#################### subroutine header end ####################

sub phylip_align {
    my $self = shift;
    my $type = shift || "dna";
    my ( $seq1, $seq2 );
    if ( $type =~ /p/i ) {
        $seq1 = $self->gaplessP1;
        $seq2 = $self->gaplessP2;
    }
    else {
        $seq1 = $self->gaplessD1;
        $seq2 = $self->gaplessD2;
    }
    my $str = "   2 " . length($seq1) . "\n";
    $str .= join( "\n", $self->name1, $seq1, $self->name2, $seq2 ) . "\n";
}

#################### subroutine header begin ####################

=head2 _check_seqs

 Usage     : $self->_check_seqs
 Purpose   : gets the protein and dna sequences for the names
 Returns   : 1 if successful for getting all four sequences, 0 otherwise
 Argument  : none
 Throws    : 0
 Comment   : used internally
           :

See Also   :

=cut

#################### subroutine header end ####################

sub _check_seqs {
    my $self = shift;

    #check to make sure that feat objects exist and are of type CDS
    unless ( $self->feat1
        && ref( $self->feat1() ) =~ /Feature/
        && $self->feat1->type->name eq "CDS" )
    {
        warn
"problem with feat1:  not defined or not a CoGeX::Feature object or Feature object is not of feature_type->name 'CDS'\n";
        return 0;
    }
    unless ( $self->feat2
        && ref( $self->feat2() ) =~ /Feature/
        && $self->feat1->type->name eq "CDS" )
    {
        warn
"problem with feat2:  not defined or not a CoGeX::Feature object or Feature object is not of feature_type->name 'CDS'\n";
        return 0;
    }
    my ($p1) = $self->feat1->protein_sequence;
    my ($p2) = $self->feat2->protein_sequence;
    my $d1   = $self->feat1->genomic_sequence;
    my $d2   = $self->feat2->genomic_sequence;
    $p1 =~ s/\s+//g;
    $p2 =~ s/\s+//g;
    if ( $p1 =~ /\*$/ ) {
        $p1 =~ s/\*$//;     #remove stop character
        $d1 =~ s/...$//;    #remove stop codon

    }
    if ( $p2 =~ /\*$/ ) {
        $p2 =~ s/\*$//;     #remove stop character
        $d2 =~ s/...$//;    #remove stop codon
    }

    #somtimes the stop aa "*" is not present while the stop codon is
    $d1 =~ s/...$// if ( length($d1) == ( 3 * length($p1) ) + 3 );
    $d2 =~ s/...$// if ( length($d2) == ( 3 * length($p2) ) + 3 );

    $self->prot1($p1)      if $p1;
    $self->prot2($p2)      if $p2;
    $self->dna1( uc($d1) ) if $d1;
    $self->dna2( uc($d2) ) if $d2;
    my ($name) = $self->feat1->names();
    $self->name1($name);
    ($name) = $self->feat2->names();
    $self->name2($name);

    return ( $self->prot1 && $self->prot2 && $self->dna1 && $self->dna2 )
      ? 1
      : 0;
}

#################### subroutine header begin ####################

=head2 _generate_gapless

 Usage     : $self->_generate_gapless
 Purpose   : called to generate the gapless alignments of the protein
             and dna sequences as determined by the aligned protein
             sequences.  This sets $self->gapless(P1|P2|D1|D3)
 Returns   : none
 Argument  : none
 Throws    : none
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub _generate_gapless {
    my $self  = shift;
    my @prot1 = split //, $self->palign1;
    my @prot2 = split //, $self->palign2;
    my ( $d1p, $d2p ) = ( 0, 0 );
    foreach my $ppos ( 0 .. $#prot1 ) {
        my $p1 = $prot1[$ppos];
        my $p2 = $prot2[$ppos];
        unless ( $p1 eq "-" || $p2 eq "-" ) {
            my $d1 = substr( $self->dna1, $d1p, 3 );
            my $d2 = substr( $self->dna2, $d2p, 3 );
            if ( $self->gaplessP1 ) {
                $self->gaplessP1( $self->gaplessP1 . $p1 );
                $self->gaplessP2( $self->gaplessP2 . $p2 );
                $self->gaplessD1( $self->gaplessD1 . $d1 );
                $self->gaplessD2( $self->gaplessD2 . $d2 );
            }
            else {
                $self->gaplessP1($p1);
                $self->gaplessP2($p2);
                $self->gaplessD1($d1);
                $self->gaplessD2($d2);
            }
        }
        #increment dna position unless in a gap
        $d1p += 3 if $p1 ne "-";
        $d2p += 3 if $p2 ne "-";
    }
}

sub _generate_DNA_alignment {
    my $self  = shift;
    my @prot1 = split //, $self->palign1;
    my @prot2 = split //, $self->palign2;
    if ( length( $self->palign1 ) ne length( $self->palign2 ) ) {
        warn
"in sub _generate_DNA_alignment.  protein sequence alignments are not of equal length!\n";
    }
    my ( $d1p, $d2p ) = ( 0, 0 );
    foreach my $ppos ( 0 .. $#prot1 ) {
        my $p1 = $prot1[$ppos];
        my $p2 = $prot2[$ppos];
        my ( $d1, $d2 );
        if ( $p1 eq "-" ) {
            $d1 = "---";
        }
        else {
            $d1 = substr( $self->dna1, $d1p, 3 );
            $d1p += 3;
        }
        if ( $p2 eq "-" ) {
            $d2 = "---";
        }
        else {
            $d2 = substr( $self->dna2, $d2p, 3 );
            $d2p += 3;
        }
        if ( $self->dalign1 ) {
            $self->dalign1( $self->dalign1 . $d1 );
            $self->dalign2( $self->dalign2 . $d2 );
        }
        else {
            $self->dalign1($d1);
            $self->dalign2($d2);
        }
    }
    return ( $self->dalign1(), $self->dalign2() );
}

1;

# The preceding line will help the module return a true value
