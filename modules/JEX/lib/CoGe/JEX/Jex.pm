package CoGe::JEX::Jex;

use v5.10;
use strict;
use warnings;

use Data::Dumper;
use Moose;
use JSON::XS;
use ZMQ::LibZMQ3;
use ZMQ::Constants qw/:all/;
use Switch;

use CoGe::JEX::Workflow;

# Attributes
has 'host' => (
    is       => 'ro',
    isa      => 'Str',
    required => 1
);

has 'port' => (
    is       => 'ro',
    isa      => 'Int',
    default  => 5151,
    required => 0
);

has 'local' => (
    is  => 'ro',
    isa => 'Bool'
);

has '_context' => (
    is       => 'ro',
    lazy     => 1,
    init_arg => undef,
    default  => sub { zmq_init(); }
);

# Public functions
sub create_workflow {
    my ( $self, %opts ) = @_;
    my $id;

    if ($opts{init} and $opts{init} == 1) {
        my $request = {
            request => 'new',
            data    => []
        };

        my $response = _send_request($self, $request);
        $id = $response->{id};

        unless ($id) {
            print STDERR "CoGe::JEX::Jex ERROR: Failed to initialize workflow ", ($opts{name} ? "'".$opts{name}."'" : '') , "\n";
            return;
        }
    }

    my $workflow = CoGe::JEX::Workflow->new(
        id      => $id,
        name    => $opts{name},
        logfile => $opts{logfile}
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
        }
    };
    
    return _send_request($self, $request);
}

sub wait_for_completion {
    my ($self, $id) = @_;
    my ($status, $wait) = (undef, 0);

    while (1) {
        $status = get_status($self, $id);

        switch ($status) {
            case /completed/i  { return 1; }
            case /notfound/i   { return 0; }
            case /failed/i     { return 0; }
            case /terminated/i { return 0; }
            case /error/i      { return 0; }
            else {
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
        }
    };

    $response = _send_request($self, $request);
    return $response->{status};
}

sub restart {
    my ($self, $id) = @_;
    my ($request, $response);

    $request = {
        request => 'restart',
        data => {
            id => $id
        }
    };

    $response = _send_request($self, $request);
    return $response->{status};
}

sub get_status {
    my ($self, $id) = @_;

    my $response = $self->get_job($id);
    return $response->{status};
}

sub get_job {
    my ($self, $id) = @_;

    my $request = {
        request => 'get_status',
        data    => {
            id => $id
        }
    };

    return _send_request($self, $request);
}

sub is_successful {
    my ($self, $response) = @_;

    return $response && lc($response->{status}) ne 'error';
}

sub get_all_workflows {
    my ($self) = @_;
    my ($request, $response, $workflows);

    $request = {
        request => 'workflows',
        data    => { }
    };

    $response = _send_request($self, $request);
    $workflows = $response->{workflows} if $response and $response->{workflows};
    $workflows //= []; #/

    return $workflows;
}

sub find_workflows {
    my ($self, $ids, $status) = @_;
    my ($request, $response, $workflows);

    my %data;
    $data{ids} = $ids if ($ids); # mdb changed 4/29/15 - rename key from "workflows" to "ids"
    $data{status} = $status if ($status);
    $request = {
        request => 'workflows',
        data    => \%data
    };

    $response = _send_request($self, $request);
    $workflows = $response->{workflows} if $response and $response->{workflows};
    $workflows //= []; #/

    return $workflows;
}

# Private functions
sub _send_request {
    my ($self, $request) = @_;
    my ($submitted, $TIMEOUT) = (time, 600);#60);#30); # mdb changed from 60 to 600 3/3/16 COGE-707

    # Set the default response as an error in case a response is not recieved.
    my $response = { status => "error" };

    eval {
        my ($socket, $msg, $json_request, $json_response);

        # mdb added 4/17/15 COGE-609 - sort json consistently
        my $json_encoder = JSON::XS->new;
        $json_encoder->canonical(1);
        
        $json_request = $json_encoder->encode($request);
        $socket = zmq_socket($self->_context, ZMQ_REQ);

        zmq_setsockopt($socket, ZMQ_LINGER, 0);
        zmq_connect($socket, _connection_string($self->host, $self->port));
        zmq_sendmsg($socket, $json_request, ZMQ_NOBLOCK);
        while ((time - $submitted) < $TIMEOUT && not defined($msg)) {
            $msg = zmq_recvmsg($socket, ZMQ_NOBLOCK);
            sleep 0.1;
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
