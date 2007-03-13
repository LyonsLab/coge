package CoGeX::Dataset;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("dataset");
__PACKAGE__->add_columns(
  "dataset_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "data_source_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "organism_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 50 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "version",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 50,
  },
  "link",
  {
    data_type => "TEXT",
    default_value => undef,
    is_nullable => 1,
    size => 65535,
  },
);

__PACKAGE__->set_primary_key("dataset_id");

__PACKAGE__->has_many("features" => "CoGeX::Feature", 'dataset_id');

__PACKAGE__->has_many("genomic_sequences" => "CoGeX::GenomicSequence", 'dataset_id');

__PACKAGE__->belongs_to("datasource" => "CoGeX::DataSource", 'data_source_id');

__PACKAGE__->belongs_to("organism" => "CoGeX::Organism", 'organism_id');

# get_sequence:
# $rs->get_genome_sequence( [chromosome, from, to] )
# $rs->get_genome_sequence() - returns the entire genome sequence *woot*
# $rs->get_genome_sequence(1, 42, 101) - returns chromosome 1 sequence 
#                                        from position 42 up to and
#                                        including position 101

sub get_genomic_sequence {
  my $self = shift;
  my %opts = @_;
  my $start = $opts{start} || $opts{begin};
  my $stop = $opts{stop} || $opts{end};
  my $chr = $opts{chr} || $opts{chromosome};
  
  my $str = "";

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
        
    $chromosome = "1" unless defined $chromosome;
    my $last = $self->last_chromosome_position($chromosome);
    $from = 1 if $from < 1;
    $to = $last if $to > $last;
    
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

sub get_genome_sequence
  {
    return shift->get_genomic_sequence(@_);
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
     my ($gs) =  $self->genomic_sequences(
					  {
					   dataset_id=>$self->dataset_id,
					   chromosome=>$chr,
					  },
					  {
					   order_by=>'stop DESC',
					   limit => 1,
					  }
					 );
     $gs->stop;
   }


1;
