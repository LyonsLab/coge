package CoGe::Exception::ItemNotFound;

use Moose;
extends 'Throwable::Error';

has message => (is => 'ro', isa => 'Str', default => 'Item not found');
has type    => (is => 'ro', isa => 'Str');
has id      => (is => 'ro', isa => 'Int');

sub as_string {
    my $self = shift;
    return $self->message . ': ' . ($self->type.' ' // '') . ($self->id // '');
}

1;