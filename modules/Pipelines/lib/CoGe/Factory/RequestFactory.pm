package CoGe::Factory::RequestFactory;

use Moose;
use CoGe::Request::Experiment;
use CoGe::Request::Genome;

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

    my $options = {
        db         => $self->db,
        jex        => $self->jex,
        user       => $self->user,
        options    => $message->{options},
        parameters => $message->{parameters}
    };

    if ($message->{type} eq "gff_export") {
        return CoGe::Request::Genome->new($options);
    }

    if ($message->{type} eq "fasta_export") {
        return CoGe::Request::Genome->new($options);
    }

    if ($message->{type} eq "experiment_export") {
        return CoGe::Request::Experiment->new($options);
    }
}

1;
