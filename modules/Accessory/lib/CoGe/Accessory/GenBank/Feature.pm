package CoGe::Accessory::GenBank::Feature;

use strict;
use Data::Dumper;
use base qw(Class::Accessor);

__PACKAGE__->mk_accessors qw(number type location qualifiers annotation blocks strand _start _stop);

sub has_qualifier {
	my $self = shift;
	my $qualifier = shift;
	if ( exists $self->qualifier->{$qualifier} ) {
		return(1);
	} else {
		return(0);
	}
}

sub get_qualifier {
	my $self = shift;
	my $qualifier = shift;
	if ( exists $self->qualifier->{$qualifier} ) {
		return($self->qualifier->{$qualifier});
	} else {
		return "?";
	}
}

sub start
  {
    my $self = shift;
    return $self->_start if defined $self->_start;
    my $start = $self->strand == 1 ? $self->_get_position() : $self->_get_position(high=>1);
    $self->_start($start);
    return $self->_start;
  }

sub stop
  {
    my $self = shift;
    return $self->_stop if defined $self->_stop;
    my $stop = $self->strand == 1 ? $self->_get_position(high=>1) : $self->_get_position();
    $self->_stop($stop);
    return $self->_stop;
  }

sub _get_position
  {
    my $self = shift;
    my %opts = @_;
    my $high = $opts{high};
    my $pos;
    foreach my $block(@{$self->blocks})
      {
	foreach my $loc (@$block)
	  {
	    $pos = $loc unless defined $pos;
	    if ($high)
	      {
		$pos = $loc if $loc > $pos;
	      }
	    else
	      {
		$pos = $loc if $loc < $pos;
	      }
	  }
      }
    return $pos;
  }

1;
