package CoGe::Factory::PipelineFactory;

use Moose;

use CoGe::Builder::GffBuilder;
use CoGe::Builder::FastaBuilder;
use CoGe::Builder::ExperimentBuilder;

has 'db' => (
    is => 'ro',
    required => 1
);

has 'conf' => (
    is => 'ro',
    required => 1
);

has 'user' => (
    is  => 'ro',
    required => 1
);

has 'jex' => (
    is => 'ro',
    required => 1
);

sub get {
    my ($self, $message) = @_;

    my $options = {
        params   => $message->{parameters},
        options  => $message->{options},
        db       => $self->db,
        jex      => $self->jex,
        user     => $self->user,
        conf     => $self->conf
    };

    my $builder;

    if ($message->{type} eq "gff_export") {
        $builder = CoGe::Builder::GffBuilder->new($options);
    } elsif ($message->{type} eq "fasta_export") {
        $builder = CoGe::Builder::FastaBuilder->new($options);
    } elsif ($message->{type} eq "experiment_export") {
        $builder = CoGe::Builder::ExperimentBuilder->new($options);
    } else {
        return;
    }

    # Construct the workflow
    $builder->build;

    # Fetch the workflow constructed
    return $builder->get;
}

1;
