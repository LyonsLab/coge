package CoGe::Exception::MissingField;

use Moose;
extends 'Throwable::Error';

has message => (is => 'ro', isa => 'Str', default => 'Missing required field');

sub as_string {
    shift->message;
}

1;