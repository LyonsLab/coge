package CoGe::Factory::RequestFactory;

use Moose;
use Data::Dumper;

use CoGe::Request::CoGeBlast;
use CoGe::Request::Empty;
use CoGe::Request::Experiment;
use CoGe::Request::NExperiments;
use CoGe::Request::ExperimentAnalysis;
use CoGe::Request::Genome;
use CoGe::Request::TwoGenomes;
use CoGe::Request::NGenomes;
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
    'blast'                 => { class => 'CoGe::Request::CoGeBlast',          authRequired => 0 },
    'export_gff'            => { class => 'CoGe::Request::Genome',             authRequired => 0 },
    'export_fasta'          => { class => 'CoGe::Request::Genome',             authRequired => 0 },
    'export_genome'         => { class => 'CoGe::Request::Genome',             authRequired => 0 },
    'export_experiment'     => { class => 'CoGe::Request::Experiment',         authRequired => 0 },
    'load_experiment'       => { class => 'CoGe::Request::Genome',             authRequired => 1 },
    'load_sra'              => { class => 'CoGe::Request::Genome',             authRequired => 1 },
    'load_batch'            => { class => 'CoGe::Request::Genome',             authRequired => 1 },
    'load_genome'           => { class => 'CoGe::Request::Empty',              authRequired => 1 },
    'load_annotation'       => { class => 'CoGe::Request::Genome',             authRequired => 1 },
    'synmap'                => { class => 'CoGe::Request::NGenomes',           authRequired => 0 },
    'synmap3d'              => { class => 'CoGe::Request::NGenomes',           authRequired => 0 },
    'pseudoassembly'        => { class => 'CoGe::Request::TwoGenomes',         authRequired => 0 },
    'dotplot_dots'          => { class => 'CoGe::Request::TwoGenomes',         authRequired => 0 },
    'analyze_snps'          => { class => 'CoGe::Request::ExperimentAnalysis', authRequired => 1 },
    'analyze_expression'    => { class => 'CoGe::Request::ExperimentAnalysis', authRequired => 1 },
    'analyze_metaplot'      => { class => 'CoGe::Request::ExperimentAnalysis', authRequired => 0 },
    'analyze_diversity'     => { class => 'CoGe::Request::ExperimentAnalysis', authRequired => 0 },
    'merge_snps'            => { class => 'CoGe::Request::NExperiments',       authRequired => 1 },
    'percent_gc_at'         => { class => 'CoGe::Request::Genome',             authRequired => 0 },
);

sub get {
    my ($self, $payload) = @_;
    unless (defined $payload && defined $payload->{type}) {
        CoGe::Exception::Generic->throw(message => 'Invalid payload', details => Dumper($payload));
    }

    my $type = $payload->{type};
    unless ($type) {
        CoGe::Exception::Generic->throw(message => "Missing job type");
    }

    my $className = $typeToClass{$type}->{class};
    unless ($className) {
        CoGe::Exception::Generic->throw(message => "Unrecognized job type: " . $payload->{type});
    }

    return $className->new(
        db           => $self->db,
        conf         => $self->conf,
        user         => $self->user,
        payload      => $payload,
        authRequired => $typeToClass{$type}->{authRequired}
    );
}

__PACKAGE__->meta->make_immutable;

1;
