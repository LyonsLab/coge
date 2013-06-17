use strict;
use warnings;
use v5.10.0;

use JSON::XS;
use ZMQ::LibZMQ3;
use ZMQ::Constants qw(ZMQ_REQ);

BEGIN {
    our($VERSION, %INTERFACES);
    
    $VERSION = 0.1;
    %INTERFACES = (
        FILE => \&_send_local,
        HTTP => \&_send_remote,
    );
}

package CoGe::Accessory::Jex;
{
    use Moose;

    # Attributes
    has 'uri' => (
        is => 'ro',
        isa  => 'Str',
        required => 1,
        handles => { uri => 'host' },
    );

    has 'local' => (
        is => 'ro',
        isa  => 'Bool',
    );

    has 'workflow' => (
        is => 'ro',
        isa => 'ArrayRef',
        writer => '_set_workflow',
    );

    has '_zmq_context' => (
        is => 'ro',
        lazy => 1,
        builder => '_build_zmq_context',
        init_arg => undef,
        handle => 'context',
    );

    sub create_workflow {
        my ($self, $name) = @_;
        my $workflow = Workflow->new(name => $name);

        return $workflow;
    }

    sub submit_workflow {
    }

    sub wait_for_compeletion {
    }

    sub terminate {
    }

    sub get_status {
    }

# Private functions
    sub _build_zmq_context {
        return zmq_init(); 
    }

    sub _send_local {
    }


    sub _send_remote {

    }
}
package CoGe::Accessory::Workflow
{
    use Moose;

    # Attributes
    has 'workflow_id' => (
        is => 'ro',
        isa  => 'Int',
        handle => 'id',
    );
    
    has 'name' => (
        is => 'ro',
        isa  => 'Str',
        handle => 'id',
    );

    has 'jobs' => (
        is => 'rw',
        isa  => 'ArrayRef',
        default => sub { [] },
    );

    # Public functions
    sub add_job {
    }
}

__PACKAGE__->meta->make_immutable();
1;
__END__
