package CoGe::Exception::Generic;

use Moose;
extends 'Throwable::Error';

has message => (is => 'ro', isa => 'Str'); # set in throw()
has details => (is => 'ro', isa => 'Str'); # set in throw()

sub as_string {
    my $self = shift;
    $self->message . ($self->details ? ': ' . $self->details : '');
}

1;
