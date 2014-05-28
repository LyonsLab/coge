package CoGe::Accessory::Jex;
use strict;
use warnings;
use v5.10;

use Moose;
use JSON::XS;
use ZMQ::LibZMQ3;
use ZMQ::Constants qw/:all/;
use CoGe::Accessory::Workflow;

# Attributes
has 'host' => (
    is       => 'ro',
    isa      => 'Str',
    required => 1,
);

has 'port' => (
    is       => 'ro',
    isa      => 'Int',
    default  => 5151,
    required => 0,
);

has 'local' => (
    is  => 'ro',
    isa => 'Bool',
);

has '_context' => (
    is       => 'ro',
    lazy     => 1,
    init_arg => undef,
    default  => sub { zmq_init(); },
);

# Public functions
sub create_workflow {
    my ( $self, %opts ) = @_;

    my $workflow = CoGe::Accessory::Workflow->new(
        id      => $opts{id},
        name    => $opts{name},
        logfile => $opts{logfile},
    );

    return $workflow;
}

sub submit_workflow {
    my ($self, $workflow) = @_;
    my ($request);

    $request = {
        request => 'schedule',
        data    => {
            id       => $workflow->id,
            name     => $workflow->name,
            logfile  => $workflow->logfile,
            priority => $workflow->priority,
            jobs     => $workflow->jobs(),
        },
    };

    return _send_request($self, $request);
}

sub wait_for_completion {
    my ($self, $id) = @_;
    my ($status, $wait) = (undef, 0);

    while (1) {
        $status = get_status($self, $id);

        given ($status) {
            when (/completed/i)  { return 1; }
            when (/notfound/i)   { return 1; }
            when (/failed/i)     { return 1; }
            when (/terminated/i) { return 2; }
            when (/error/i)      { return 3; }
            default {
                sleep $wait;
                $wait = $wait + 0.25;
            }
        }
    }
}

sub terminate {
    my ($self, $id) = @_;
    my ($request, $response);

    $request = {
        request => 'cancel',
        data => {
            id => $id
        },
    };

    $response = _send_request($self, $request);
    return $response->{status};
}

sub get_status {
    my ($self, $id) = @_;
    my ($request, $response);

    $request = {
        request => 'get_status',
        data    => {
            id => $id
        },
    };

    $response = _send_request($self, $request);
    return $response->{status};
}

sub get_all_workflows {
    my ($self) = @_;
    my ($request, $response, $workflows);

    $request = {
        request => 'get_workflows',
        data    => {},
    };

    $response = _send_request($self, $request);
    $workflows = $response->{workflows} if $response and $response->{workflows};
    $workflows //= [];

    return $workflows;
}

# Private functions
sub _send_request {
    my ($self, $request) = @_;
    my ($waited, $TIMEOUT) = (0, 30);

    # Set the default response as an error in case a response is not recieved.
    my $response = { status => "error" };

    eval {
        my ($socket, $msg, $json_request, $json_response);

        $json_request = encode_json($request);
        $socket = zmq_socket($self->_context, ZMQ_REQ);

        zmq_connect($socket, _connection_string($self->host, $self->port));
        zmq_sendmsg($socket, $json_request, ZMQ_NOBLOCK);

        while ($TIMEOUT > $waited && not defined($msg)) {
            $msg = zmq_recvmsg($socket, ZMQ_NOBLOCK);
            $waited++;
            sleep 1;
        }

        zmq_close($socket);

        $json_response = zmq_msg_data($msg) if $msg;
        $response = decode_json($json_response) if $json_response;
    };

    # Make sure an error did not occur
    die $@ if $@;

    return $response;
}

sub _connection_string {
    my ($host, $port) = @_;
    return "tcp://$host:$port";
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
__END__

