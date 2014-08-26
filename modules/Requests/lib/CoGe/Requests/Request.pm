package CoGe::Requests::Request;

use Moose::Role;

has 'options' => (
    is        => 'ro',
    isa       => 'HashRef',
    required  => 1
);

has 'parameters' => (
    is        => 'ro',
    isa       => 'HashRef',
    required  => 1
);

requires qw(is_valid has_access execute);

1;
