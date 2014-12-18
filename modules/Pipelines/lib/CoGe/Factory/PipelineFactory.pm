package CoGe::Factory::PipelineFactory;

use Moose;

use CoGe::Builder::Export::Gff;
use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Experiment;

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

    # Select pipeline builder
    my $builder;
    if ($message->{type} eq "export_gff") {
        $builder = CoGe::Builder::Export::Gff->new($options);
    } 
    elsif ($message->{type} eq "export_fasta") {
        $builder = CoGe::Builder::Export::Fasta->new($options);
    } 
    elsif ($message->{type} eq "export_experiment") {
        $builder = CoGe::Builder::Export::Experiment->new($options);
    } 
    else {
        return;
    }

    # Construct the workflow
    $builder->build;

    # Fetch the workflow constructed
    return $builder->get;
}

1;
