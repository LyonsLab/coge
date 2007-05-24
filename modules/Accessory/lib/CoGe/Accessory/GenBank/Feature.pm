package CoGe::Accessory::GenBank::Feature;

use strict;
use Data::Dumper;
use base qw(Class::Accessor);

__PACKAGE__->mk_accessors qw(number type location qualifiers annotation blocks);

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

1;
