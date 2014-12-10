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

has 'user'  => (
    is        => 'ro',
    required  => 1
);

has 'db' => (
    is        => 'ro',
    required  => 1
);

has 'jex' => (
    is        => 'ro',
    required  => 1
);

our $NOT_FOUND = "the job could be found";

sub execute {
    my ($self, $workflow) = @_;

    # Check for workflow
    return {
        success => JSON::false,
        error => { Invalid => $NOT_FOUND }
    } unless $workflow;

    my $resp = $self->jex->submit_workflow($workflow);
    my $success = $self->jex->is_successful($resp);

    return {
        job_id => $resp->{id},
        success => $success ? JSON::true : JSON::false
    };
}

1;
