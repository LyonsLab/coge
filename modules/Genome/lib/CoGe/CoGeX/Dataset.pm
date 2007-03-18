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



1;
