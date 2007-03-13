package CoGe::CoGeX::Dataset;


use strict;
use warnings;

use base 'CoGeX::Dataset';

 sub data_source
   {
     my $self = shift;
     return $self->datasource();
   }

 sub source
   {
     my $self = shift;
     return $self->datasource();
   }

 sub feature
   {
     my $self = shift;
     return $self->features();
     }

 sub feats
   {
     my $self = shift;
     return $self->features();
     }

 sub genomic_sequences
   {
     my $self = shift;
     return $self->genome_sequences();
   }

 sub genomic_seqs
   {
     my $self = shift;
     return $self->genome_sequences();
   }

 sub seqs
   {
     my $self = shift;
     return $self->genome_sequence();
   }

 sub desc
   {
     my $self = shift;
     return $self->description();
   }

 sub org
   {
     my $self = shift;
     return $self->organism();
   }

 sub species
   {
     my $self = shift;
     return $self->organism();
   }

 sub url
   {
     my $self = shift;
     return $self->link();
   }

 sub v
   {
     my $self = shift;
     return $self->version();
   }

 sub ver
   {
     my $self = shift;
     return $self->version();
   }


################################################## subroutine header start ##

=head2 last_chromsome_position

 Usage     : my $last = $genome_seq_obj->last_chromosome_position($chr);
 Purpose   : gets the last genomic sequence position for a dataset given a chromosome
 Returns   : an integer that refers to the last position in the genomic sequence refered
             to by a dataset given a chromosome
 Argument  : string => chromsome for which the last position is sought
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


 sub last_chromosome_position
   {
     my $self = shift;
     my $chr = shift;
     return $self->genomic_sequences(
				     {
				      dataset_id=>$self->id,
				      chromosome=>$chr,
				     },
				     {
				      order_by=>'stop DESC',
				      limit => 1,
				     }
				    );
				   
   }

 # get_sequence:
 # $rs->get_genome_sequence( [chromosome, from, to] )
 # $rs->get_genome_sequence() - returns the entire genome sequence *woot*
 # $rs->get_genome_sequence(1, 42, 101) - returns chromosome 1 sequence 
 #                                        from position 42 up to and
 #                                        including position 101

 sub get_genome_sequence {
   my $self = shift;
   my %opts = @_;
   my $start = $opts{start} || $opts{begin};
   my $stop = $opts{stop} || $opts{end};
   my $chr = $opts{chr} || $opts{chromosome};
   my $str = "";
   use Data::Dumper;
   print STDERR Dumper \%opts;
   if ( @_ > 1 ) {
     my($chromosome, $from, $to);
     if (defined $start && defined $stop && defined $chr)
       {
 	$chromosome = $chr;
 	$from = $start;
 	$to = $stop;
       }
     else
       {
 	($chromosome, $from, $to) = @_;
       }
     print STDERR join ("\t", $self->name, $chromosome, $from, $to),"\n";

     $chromosome = "1" unless defined $chromosome;
     my $last = $self->last_chromosome_position($chromosome);
     $from = 1 if $from < 1;
     $to = $last if $to > $last;
     print STDERR join ("\t", $self->name, $chromosome, $from, $to),"\n";
     # make sure two numbers were sent in
     return undef unless ($from =~ /\A\d+\z/ and  $to =~ /\A\d+\z/);
     return undef unless $to > $from;
     my $g1 = $self->genome_sequences(
                       {
                         'chromosome' => $chromosome,
                         'start'      => { '<=' => $from },
                         'stop'       => { '>=' => $from }
                       })->first();
     my $g2 = $self->genome_sequences(
                       {
                         'chromosome' => $chromosome ,
                         'start'      => { '<=' => $to },
                         'stop'       => { '>=' => $to }
                       })->first();
 #    use Data::Dumper;
 #    print STDERR Dumper ($g1,$g2);
     if ( $g2->start == $g1->start )         { # hit is within one row
       $str = $g1->sequence_data;
     } elsif ( $g2->start == $g1->stop + 1 ) { # consecutive rows
       $str = $g1->sequence_data . $g2->sequence_data;
     } elsif ( $g2->start >= $g1->stop + ($g1->stop - $g1->start) ) {
       $str = $g1->sequence_data();            # start with first row
       my $inneriter = $self->genome_sequences(
                                   {
                                     'chromosome' => $chromosome ,
                                     'start' => { '>=' => $g1->stop + 1 },
                                     'stop ' => { '<'  => $g2->start }
                                   });

        
       while( my $innerseq = $inneriter->next() ) {
         $str .= $innerseq->sequence_data();
       }
       $str .= $g2->sequence_data(); # get last row
     }
     $str = $self->trim_sequence( $str, $g1->start, $g2->stop, $from, $to );

   } elsif ( @_ == 1 ) {    # get a whole chromosome
     my $chromosome = shift;
     my $allseqs = $self->genome_sequences( { 'chromosome' => $chromosome } );
     while ( my $g = $allseqs->next ) {
       $str .= $g->sequence_data;
     }
   } else {                 # entire sequence
     my $allseqs = $self->genome_sequences();
     while ( my $g = $allseqs->next ) {
       $str .= $g->sequence_data;
     }
   }
   return $str;
 }

#         .         .         .         .
#1234567890123456789012345678901234567890
#                                       CCACAACCAGCTGACTAGGTA
#ACGACGCAGCTATGGCCTCCCCGCCCACCAGGCCGCCAGCCACAACCAGC
#         CTATGGCCTC
sub trim_sequence {
  my $self = shift;
  my( $seq, $seqstart, $seqend, $newstart, $newend ) = @_;
  # substr( string, offset, length ) will leave off "length"
  # many characters from the end of the string if length is negative
  $seq = substr($seq, $newstart - $seqstart, - ($seqend - $newend) );
  return($seq);
}

1;
