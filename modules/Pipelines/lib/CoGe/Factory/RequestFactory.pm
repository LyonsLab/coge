package CoGe::Factory::RequestFactory;

use Moose;
use CoGe::Requests::GffRequest;

has 'user'    => (
    is        => 'ro',
    required  => 1
);

has 'db'      => (
    is        => 'ro',
    required  => 1
);

has 'jex'     => (
    is        => 'ro',
    required  => 1
);

sub get {
    my ($self, $message) = @_;

    if ($message->{type} eq "gff_export") {
        return CoGe::Requests::GffRequest->new(
            db         => $self->db,
            jex        => $self->jex,
            user       => $self->user,
            options    => $message->{options},
            parameters => $message->{parameters}
        );
    }
}

1;
