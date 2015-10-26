package CoGe::Builder::Buildable;

use Moose::Role;

requires qw(build);

has 'workflow' => (
    is => 'rw',
    isa  => 'CoGe::Accessory::Workflow'
);

# mdb removed 10/12/15 - no longer needed
#has 'jex' => (
#    is => 'ro'
#);

has 'staging_dir' => (
    is => 'rw'
);

has 'result_dir' => (
    is => 'rw'
);

has 'params' => (
    is => 'ro',
    required => 1
);

has 'requester' => (
    is => 'ro'
);

has 'site_url' => ( # mdb added 10/12/15
    is => 'rw'
);

has 'page' => ( # mdb added 10/12/15
    is => 'rw'
);

has 'db' => (
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

has 'outputs' => (
    is => 'rw',
    default => sub { {} } # initialize to empty hash
);

# mdb removed 10/12/15 - doesn't make sense
#sub get {
#    my $self = shift;
#    return $self->workflow;
#}

1;
