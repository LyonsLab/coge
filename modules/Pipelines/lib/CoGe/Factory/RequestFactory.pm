package CoGe::Factory::RequestFactory;

use Moose;
use CoGe::Request::Experiment;
use CoGe::Request::ExperimentAnalysis;
use CoGe::Request::Genome;
use CoGe::Request::Empty;

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
    unless (defined $message && defined $message->{type}) {
        print STDERR "RequestFactory: error: invalid message\n";
        return;
    }

    my $options = {
        db         => $self->db,
        jex        => $self->jex,
        user       => $self->user,
        parameters => $message->{parameters}
    };

    my $type = $message->{type};
    if ($type eq "export_gff" ||
        $type eq "export_fasta" ||
        $type eq "export_genome" ||
        $type eq "load_experiment" ||
        $type eq "load_batch" ||
        $type eq "load_annotation")
    {
        return CoGe::Request::Genome->new($options);
    }
    elsif ($type eq "export_experiment") {
        return CoGe::Request::Experiment->new($options);
    }
    elsif ($type eq "analyze_snps") {
        return CoGe::Request::ExperimentAnalysis->new($options);
    }
    elsif ($type eq "load_genome") {
        return CoGe::Request::Empty->new($options);
    }
    else {
        print STDERR "RequestFactory: error: unrecognized job type '", $type, "'\n";
        return;
    }
}

1;
