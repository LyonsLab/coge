package CoGe::Accessory::Workflow;
use strict;
use warnings;
use 5.10.0;

use Moose;

# Attributes
has 'workflow_id' => (
    is => 'ro',
);

has 'name' => (
    is => 'ro',
    isa  => 'Str',
    required => 1,
    default => "",
);

has 'jobs' => (
    is => 'rw',
    isa  => 'ArrayRef',
    default => sub { [] },
    required => 0,
);

# Public functions
sub add_job {
    my ($self, %opts) = @_;
    my $cmd = $opts{cmd};
    my $script = "" unless defined($opts{script});
    my $args = $opts{args};
    my $inputs = $opts{inputs};
    my $outputs = $opts{outputs};

    my $size = $self->jobs;

    push(@{$self->jobs}, {
        cmdstring => {cmd => $cmd, script => $script, args => $args,},
        inputs => $inputs,
        outputs => $outputs,
    });

    return scalar($self->jobs()) > $size;
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
__END__
