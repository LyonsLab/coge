package CoGe::Request::Request;

use Moose;
use Data::Dumper;

has 'parameters' => (
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

has 'jex' => (
    is        => 'ro',
    required  => 1
);

sub execute {
    my ($self, $pipeline) = @_;

    # Check for workflow
    my $workflow = $pipeline->workflow;
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

    my $response = {
        id => $resp->{id},
        success => JSON::true        
    };
    $response->{site_url} = $pipeline->site_url if ($pipeline->site_url);
    return $response;
}

__PACKAGE__->meta->make_immutable;

1;
