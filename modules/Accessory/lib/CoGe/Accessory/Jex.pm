package CoGe::Accessory::Jex;
use strict;
use warnings;
use 5.10.0;

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
        id    => $opts{id},
        name  => $opts{name},
        logfile => $opts{logfile},
    );

    return $workflow;
}

sub submit_workflow {
    my ( $self, $workflow ) = @_;
    my ( $socket, $msg, $jobs );
    $jobs = $workflow->jobs();

    my $request = encode_json(
        {
            request => 'schedule',
            data    => {
                id      => $workflow->id,
                name    => $workflow->name,
                logfile => $workflow->logfile,
                jobs    => $jobs,
            },
        }
    );

    $socket = zmq_socket( $self->_context, ZMQ_REQ );
    zmq_setsockopt($socket, ZMQ_LINGER, 0);
    zmq_connect( $socket, _connection_string( $self->host, $self->port ) );

    zmq_sendmsg( $socket, $request, ZMQ_NOBLOCK);

    my $count = 0;

    while (5 > $count && not defined($msg)) {
        $msg = zmq_recvmsg($socket, ZMQ_NOBLOCK);
        $count++;
        sleep 1;
    }

    zmq_close($socket);
    zmq_term($self->_context);

    my $result = zmq_msg_data($msg) if $msg;
    $result //= '{"error" : "ERROR"}';

    return $result;
}

sub wait_for_completion {
    my ( $self, $id ) = @_;
    my ( $status, $wait );
    $wait = 0;

    while (1) {
        $status = get_status( $self, $id );

        given ($status) {
            when ("Completed")  { return 1; }
            when ("NotFound")   { return 1; }
            when ("Failed")     { return 1; }
            when ("Terminated") { return 2; }
            when ("Error")      { return 3; }
            default {
                sleep $wait;
                $wait = $wait + 0.25;
            }
        }
    }
}

sub terminate {
    my ( $self, $id ) = @_;
    my ( $socket, $resp, $cmd, $msg );

    $cmd = encode_json(
        {   request => 'cancel',
            data => {
                id => $id
            },
        }
    );
    $socket = zmq_socket( $self->_context, ZMQ_REQ );

    zmq_connect( $socket, _connection_string( $self->host, $self->port ) );
    zmq_send( $socket, $cmd );

    $resp = zmq_recvmsg($socket);
    $msg  = decode_json( zmq_msg_data($resp) );

    return $msg;
}

sub get_status {
    my ($self, $id ) = @_;
    my ( $socket, $reply, $msg );

    $socket = zmq_socket( $self->_context, ZMQ_REQ );
    zmq_connect( $socket, _connection_string( $self->host, $self->port ) );

    my $request = encode_json(
        {
            request => 'get_status',
            data    => {
                id => $id
            },
        }
    );

    $msg = zmq_send( $socket, $request );
    $msg = zmq_recvmsg($socket);

    my $data = zmq_msg_data($msg);
    my $res  = decode_json($data);

    return ${$res}{'status'};
}

# Private functions
sub _build_zmq_context {
    return zmq_init();
}

sub _connection_string {
    my ( $host, $port ) = @_;
    return "tcp://$host:$port";
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
__END__

