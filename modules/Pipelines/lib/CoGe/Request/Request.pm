package CoGe::Request::Request;

use Moose::Role;

use Data::Dumper;

has 'options' => (
    is        => 'ro',
    #isa       => 'HashRef',
    #required  => 1
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
        error => { Error => "failed to build workflow" }
    } unless $workflow;

    my $resp = $self->jex->submit_workflow($workflow);
    my $success = $self->jex->is_successful($resp);
    
    unless ($success) {
        print STDERR 'JEX response: ', Dumper $resp, "\n";
        return {
            success => JSON::false,
            error => { JEX => 'failed to submit workflow' }
            #TODO return $resp error message from JEX
        }
    }

    return {
        id => $resp->{id},
        success => JSON::true
    };
}

1;
