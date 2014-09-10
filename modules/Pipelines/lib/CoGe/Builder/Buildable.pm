package CoGe::Builder::Buildable;

use Moose::Role;

requires qw(build);

has 'workflow' => (
    is => 'rw',
    isa  => 'CoGe::Accessory::Workflow'
);

has 'options' => (
    is => 'ro',
    required => 1
);

has 'params' => (
    is => 'ro',
    required => 1
);

has 'db' => (
    is => 'ro',
    required => 1
);

has 'jex' => (
    is => 'ro',
    required => 1
);

has 'user' => (
    is => 'ro',
    required => 1
);

has 'conf' => (
    is => 'ro',
    required => 1
);

sub get {
    my $self = shift;
    return $self->workflow;
}

1;
