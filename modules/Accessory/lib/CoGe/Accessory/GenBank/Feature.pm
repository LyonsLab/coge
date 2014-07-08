package CoGe::Accessory::GenBank::Feature;

use strict;
use Data::Dumper;
use base qw(Class::Accessor);

__PACKAGE__->mk_accessors('number', 'type', 'location', 'qualifiers',
	'annotation', 'blocks', 'strand', '_start', '_stop');

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
    my $start = $self->_get_position();
    $self->_start($start);
    return $self->_start;
  }

sub stop
  {
    my $self = shift;
    return $self->_stop if defined $self->_stop;
    my $stop = $self->_get_position(high=>1);
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

=head1 NAME

Feature

=head1 SYNOPSIS

use Feature

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
