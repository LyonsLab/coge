<<<<<<< HEAD
package CoGe::Exception::Generic;
=======
package CoGe::Exception::Custom;
>>>>>>> 64ad5e911b0cc224edca76011b04410d5f95858d

use Moose;
extends 'Throwable::Error';

has message => (is => 'ro', isa => 'Str'); # set in throw()
<<<<<<< HEAD
has details => (is => 'ro', isa => 'Str'); # set in throw()

sub as_string {
    my $self = shift;
    return $self->message . ($self->details ? ': ' . $self->details : '');
}

1;
=======

sub as_string {
    shift->message;
}

1;
>>>>>>> 64ad5e911b0cc224edca76011b04410d5f95858d
