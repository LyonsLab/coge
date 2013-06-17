package CoGe::Accessory::Jex;
use strict;
use warnings;
use 5.10.0;

use Moose;
use JSON::XS;
use ZMQ::LibZMQ2;
use ZMQ::Constants qw(ZMQ_REQ);
use CoGe::Accessory::Workflow;


# Attributes
has 'host' => (
    is => 'ro',
    isa  => 'Str',
    required => 1,
);

has 'port' => (
    is => 'ro',
    isa => 'Int',
    default => 5151,
    required => 0,
);

has 'local' => (
    is => 'ro',
    isa  => 'Bool',
);

has '_context' => (
    is => 'ro',
    lazy => 1,
    init_arg => undef,
    default => sub { zmq_init(); },
);

# Public functions
sub create_workflow {
    my ($self, $name) = @_;
    my $workflow = CoGe::Accessory::Workflow->new(name => $name);

    return $workflow;
}

sub submit_workflow {
    my ($self, $workflow) = @_;
    my ($socket, $msg, $jobs);
    $jobs = $workflow->jobs();

    my $request = encode_json({
        request => 'schedule',
        data => {
                    'name' => $workflow->name,
                    'jobs' => $jobs,
                },
    });

    $socket = zmq_socket($self->_context, ZMQ_REQ);
    zmq_connect($socket, _connection_string($self->host, $self->port));
    zmq_send($socket, $request);
    $msg = zmq_recv($socket);

    return zmq_msg_data($msg);
}

sub wait_for_completion {
    my ($self, $id) = @_;
    my ($status, $wait);
    $wait = 1;
    
    while (1) {
        $status = get_status($self, $id);

        say $status;
        if ($status eq "TERMINATED") {
            return -1;
        }

        if ($status eq "COMPLETED") {
            return 1;
        }
        
        if ($status eq "CANCELLED") {
            return -1;
        }
        
        if ($status eq "NOT_FOUND") {
            return -1;
        }

        sleep $wait;
        $wait = $wait * 2;
    }
}

sub terminate {
    my ($self, $id) = @_;
    my ($socket, $resp, $cmd, $msg);

    $cmd = encode_json { cmd => 'TERM', id => $id};
    $socket = zmq_socket($self->_context, ZMQ_REQ);

    zmq_connect($socket, _connection_string($self->host, $self->port));
    zmq_send($socket, $cmd);

    $resp = zmq_recv($socket);
    $msg = decode_json(zmq_msg_data($resp));

    return $msg;
}

sub get_status {
    my ($self, $id) = @_;
    my ($socket, $reply, $msg);

    $socket = zmq_socket($self->_context, ZMQ_REQ);
    zmq_connect($socket, _connection_string($self->host, $self->port));

    my $request = encode_json({
        request => 'get_status',
        data => $id,
    });
    
    $msg = zmq_send($socket, $request);
    $msg = zmq_recv($socket);

    my $data = zmq_msg_data($msg);
    my $res = decode_json($data);

    return ${$res}{'status'};
}

# Private functions
sub _build_zmq_context {
    return zmq_init();
}

sub _connection_string {
    my ($host, $port) = @_;
    return "tcp://$host:$port";
}

sub _send_local {
}


sub _send_remote {

}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
__END__

