package CoGe::Factory::PipelineFactory;

use Moose;

use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Gff;
use CoGe::Builder::Export::Experiment;
use CoGe::Builder::Load::Experiment;
use CoGe::Builder::Load::BatchExperiment;
use CoGe::Builder::Load::Genome;
use CoGe::Builder::Load::Annotation;
use CoGe::Builder::SNP::IdentifySNPs;

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

    my $request = {
        params    => $message->{parameters},
        requester => $message->{requester},
        db        => $self->db,
        jex       => $self->jex,
        user      => $self->user,
        conf      => $self->conf
    };

    # Select pipeline builder
    my $builder;
    if ($message->{type} eq "export_gff") {
        $builder = CoGe::Builder::Export::Gff->new($request);
    } 
    elsif ($message->{type} eq "export_fasta") {
        $builder = CoGe::Builder::Export::Fasta->new($request);
    } 
    elsif ($message->{type} eq "export_experiment") {
        $builder = CoGe::Builder::Export::Experiment->new($request);
    }
    elsif ($message->{type} eq "load_experiment") {
        $builder = CoGe::Builder::Load::Experiment->new($request);
    }
    elsif ($message->{type} eq "load_batch") {
        $builder = CoGe::Builder::Load::BatchExperiment->new($request);
    }
    elsif ($message->{type} eq "load_genome") {
        $builder = CoGe::Builder::Load::Genome->new($request);
    }
    elsif ($message->{type} eq "load_annotation") {
        $builder = CoGe::Builder::Load::Annotation->new($request);
    }
    elsif ($message->{type} eq "analyze_snps") {
        $builder = CoGe::Builder::SNP::IdentifySNPs->new($request);
    }
    else {
        print STDERR "PipelineFactory::get unknown type\n";
        return;
    }

    # Construct the workflow
    my $rc = $builder->build;
    unless ($rc) {
        print STDERR "PipelineFactory::get build failed\n";
        return;
    }

    # Fetch the workflow constructed
    return $builder->get;
}

1;
