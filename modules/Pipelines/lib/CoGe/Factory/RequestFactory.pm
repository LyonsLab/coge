package CoGe::Factory::RequestFactory;

use Moose;

use CoGe::Request::CoGeBlast;
use CoGe::Request::Empty;
use CoGe::Request::Experiment;
use CoGe::Request::ExperimentAnalysis;
use CoGe::Request::Genome;
use CoGe::Request::SynMap;
use CoGe::Request::TwoGenomes;

use Data::Dumper;
use Switch;

has 'user'    => (
    is        => 'ro',
    required  => 1
);

has 'db'      => (
    is        => 'ro',
    required  => 1
);

has 'conf' => (
    is => 'ro',
    required => 1
);

has 'jex'     => (
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
    'analyze_snps'          => 'CoGe::Request::ExperimentAnalysis',
    'synmap'                => 'CoGe::Request::SynMap',
    'synmap3d'              => 'CoGe::Request::SynMap',
    'dotplot_dots'          => 'CoGe::Request::TwoGenomes',
    'analyze_expression'    => 'CoGe::Request::ExperimentAnalysis',
    'analyze_metaplot'      => 'CoGe::Request::ExperimentAnalysis',
    'analyze_diversity'     => 'CoGe::Request::ExperimentAnalysis',
    'nuccounter'            => 'CoGe::Request::Genome'
);

sub get {
    my ($self, $message) = @_;
    unless (defined $message && defined $message->{type}) {
        warn 'RequestFactory: error: invalid message';
        warn Dumper $message;
        return;
    }

    my $className = $typeToClass{$message->{type}};
    unless ($className) {
        warn "RequestFactory: error: unrecognized job type '", $message->{type};
        return;
    }

    return $className->new(
        db         => $self->db,
        conf       => $self->conf,
        jex        => $self->jex,
        user       => $self->user,
        parameters => $message->{parameters}
    };
}

1;
