package CoGe::Request::Request;

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

sub execute {
    my ($self, $workflow) = @_;

    # Check for workflow
    return {
        success => JSON::false,
        error => { Invalid => "job not found" }
    } unless $workflow;

    my $resp = $self->jex->submit_workflow($workflow);
    my $success = $self->jex->is_successful($resp);

    return {
        id => $resp->{id},
        success => $success ? JSON::true : JSON::false
    };
}

1;
