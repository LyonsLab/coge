package CoGe::Request::Request;

use Moose;
use Data::Dumper;

has 'payload' => (
    is        => 'ro',
    isa       => 'HashRef',
    required  => 1
);

has 'user'  => (
    is        => 'ro',
    isa       => 'Maybe[CoGeX::Result::User]', # "Maybe" is to prevent error on undef
    required  => 1
);

has 'db' => (
    is        => 'ro',
    isa       => 'CoGeX',
    required  => 1
);

has 'conf' => (
    is       => 'ro',
    required => 1
);

sub type       { shift->payload->{type} }
sub parameters { shift->payload->{parameters} }
sub requester  { shift->payload->{requester} }

__PACKAGE__->meta->make_immutable;

__PACKAGE__->meta->make_immutable;

1;
