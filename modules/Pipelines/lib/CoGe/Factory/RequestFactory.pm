package CoGe::Factory::RequestFactory;

use Moose;
use CoGe::Requests::ExperimentRequest;
use CoGe::Requests::GenomeRequest;
use Data::Dumper;

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
    return unless $message;
    #print STDERR Dumper $message, "\n";

    my $options = {
        db         => $self->db,
        jex        => $self->jex,
        user       => $self->user,
        options    => $message->{options},
        parameters => $message->{parameters}
    };

    if ($message->{type} eq "export_gff") {
        return CoGe::Requests::GenomeRequest->new($options);
    }
    if ($message->{type} eq "export_fasta") {
        return CoGe::Requests::GenomeRequest->new($options);
    }
    if ($message->{type} eq "export_experiment") {
        return CoGe::Requests::ExperimentRequest->new($options);
    }
    if ($message->{type} eq "export_genome") {
        return CoGe::Requests::GenomeRequest->new($options);
    }
}

1;
