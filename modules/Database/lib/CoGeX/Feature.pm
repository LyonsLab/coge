package CoGeX::Feature;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("feature");
__PACKAGE__->add_columns(
  "feature_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "feature_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "dataset_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
);

__PACKAGE__->set_primary_key("feature_id");

# feature has many feature_names
__PACKAGE__->has_many( 'feature_names' => "CoGeX::FeatureName", 'feature_id');

# feature has many annotations
__PACKAGE__->has_many( 'annotations' => "CoGeX::Annotation", 'feature_id');

# feature has many locations
__PACKAGE__->has_many( 'locations' => "CoGeX::Location", 'feature_id');

# feature has many sequences - note this is non-genomic sequence (see
# genomic_sequence() and sequence
__PACKAGE__->has_many( 'sequences' => "CoGeX::Sequence", 'feature_id');

# feature_type has many features
__PACKAGE__->belongs_to("feature_type" => "CoGeX::FeatureType", 'feature_type_id');

# dataset has many features
__PACKAGE__->belongs_to("dataset" => "CoGeX::Dataset", 'dataset_id');

sub esearch : ResultSet {
    my $self = shift;
    $_[1]{'join'} = ['feature_type','feature_names','annotations','locations'],
    $_[1]{'prefetch'} = ['feature_type','feature_names'],
    return $self->search(
         @_ 
    );

}

sub protein_sequence {
  my $self = shift;
  my $siter = $self->sequences();
  my @sequence_objects;
  while ( my $seq = $siter->next() ) {
    push @sequence_objects, $seq;
  }

  if (@sequence_objects == 1) {
    return $sequence_objects[0]->sequence_data();
  } elsif ( @sequence_objects > 1 )  {
    return \@sequence_objects;
  } else {
    return undef;
  }
}

sub genome_sequence {
  my $self = shift;
  my $dataset = $self->dataset();
  my @sequences;
  my $strand = "+";
  my $lociter = $self->locations();
  while ( my $loc = $lociter->next() ) {
    my $fseq = $dataset->get_genome_sequence(
                                             $loc->chromosome(),
                                             $loc->start,
                                             $loc->stop );
    if ( $loc->strand eq "-" or $loc->strand eq "-1" ) {
      push @sequences, $self->revcomp($fseq);
    } else {
      push @sequences, $fseq;
    }
  }
  return wantarray ? @sequences : join( "", @sequences );
}

sub reverse_comp {
  my $self = shift;
  my $seq = shift;
  my $revseq = uc(reverse $seq);
  $revseq =~ tr/ATCG/TAGC/;
  return $revseq;
}

sub revcomp {
  my $self = shift;
  return $self->reverse_comp( @_ );
}

1;
