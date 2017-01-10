package CoGe::Exception::Custom;

use Moose;
extends 'Throwable::Error';

has message => (is => 'ro', isa => 'Str'); # set in throw()

sub as_string {
    shift->message;
}

1;