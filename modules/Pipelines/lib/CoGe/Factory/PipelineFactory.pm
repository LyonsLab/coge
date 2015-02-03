package CoGe::Factory::PipelineFactory;

use Moose;

use CoGe::Builder::GffBuilder;
use CoGe::Builder::FastaBuilder;
use CoGe::Builder::ExperimentBuilder;
use CoGe::Builder::GenomeBuilder;

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
    if ($message->{type} eq "export_gff") {
        $builder = CoGe::Builder::GffBuilder->new($options);
    } 
    elsif ($message->{type} eq "export_fasta") {
        $builder = CoGe::Builder::FastaBuilder->new($options);
    } 
    elsif ($message->{type} eq "export_experiment") {
        $builder = CoGe::Builder::ExperimentBuilder->new($options);
    }
    elsif ($message->{type} eq "export_genome") {
        $builder = CoGe::Builder::GenomeBuilder->new($options);
    } 
    else {
        print STDERR "CoGe::Factory::PipelineFactory: unknown job type '", ($message->{type} ? $message->{type} : ''), "'\n";
        return; # error: unknown job type
    }

    # Construct the workflow
    unless ($builder->build) {
        print STDERR "CoGe::Factory::PipelineFactory: build failed\n";
        return; # error: build failed
    }

    # Fetch the workflow constructed
    return $builder->get;
}

1;
