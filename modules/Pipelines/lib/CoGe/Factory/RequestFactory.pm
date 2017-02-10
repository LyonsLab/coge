package CoGe::Factory::RequestFactory;

use Moose;
use Data::Dumper;

use CoGe::Request::CoGeBlast;
use CoGe::Request::Empty;
use CoGe::Request::Experiment;
use CoGe::Request::ExperimentAnalysis;
use CoGe::Request::Genome;
use CoGe::Request::TwoGenomes;
use CoGe::Request::NGenomes; # TODO: replace Genome and TwoGenomes with this
use CoGe::Exception::Generic;

has 'user'    => (
    is        => 'ro',
    required  => 1
);

has 'db'      => (
    is        => 'ro',
    required  => 1
);

has 'conf'    => (
    is        => 'ro',
    required  => 1
);

my %typeToClass = (
    'blast'                 => 'CoGe::Request::CoGeBlast',
    'export_gff'            => 'CoGe::Request::Genome',
    'export_fasta'          => 'CoGe::Request::Genome',
    'export_genome'         => 'CoGe::Request::Genome',
    'export_experiment'     => 'CoGe::Request::Experiment',
    'load_experiment'       => 'CoGe::Request::Genome',
    'load_sra'              => 'CoGe::Request::Genome',
    'load_batch'            => 'CoGe::Request::Genome',
    'load_genome'           => 'CoGe::Request::Empty',
    'load_annotation'       => 'CoGe::Request::Genome',
    'synmap'                => 'CoGe::Request::NGenomes',
    'synmap3d'              => 'CoGe::Request::NGenomes',
    'pseudoassembly'        => 'CoGe::Request::TwoGenomes',
    'dotplot_dots'          => 'CoGe::Request::TwoGenomes',
    'analyze_snps'          => 'CoGe::Request::ExperimentAnalysis',
    'analyze_expression'    => 'CoGe::Request::ExperimentAnalysis',
    'analyze_metaplot'      => 'CoGe::Request::ExperimentAnalysis',
    'analyze_diversity'     => 'CoGe::Request::ExperimentAnalysis',
    'analyze_nucleotides'   => 'CoGe::Request::Genome'
);

sub get {
    my ($self, $payload) = @_;
    unless (defined $payload && defined $payload->{type}) {
        CoGe::Exception::Generic->throw(message => 'Invalid payload', details => Dumper($payload));
    }

    my $className = $typeToClass{$payload->{type}};
    unless ($className) {
        CoGe::Exception::Generic->throw(message => "Unrecognized job type: " . $payload->{type});
    }

    return $className->new(
        db      => $self->db,
        conf    => $self->conf,
        user    => $self->user,
        payload => $payload
    );
}

__PACKAGE__->meta->make_immutable;

1;
