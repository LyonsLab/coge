package CoGe::CoGeX::Location;

use strict;
use warnings;

use base 'CoGeX::Location';

sub begin
  {
    my $self = shift;
    return $self->start(@_);
  }

sub end
  {
    my $self = shift;
    return $self->stop(@_);
  }

sub chr
  {
    my $self = shift;
    return $self->chromosome(@_);
  }

sub feat
  {
    my $self = shift;
    return $self->feature();
  }


1;

