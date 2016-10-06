package CoGe::Builder::Buildable;

use Moose::Role;

requires qw(build);

has 'workflow' => (
    is => 'rw',
    isa  => 'CoGe::JEX::Workflow'
);

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

1;
