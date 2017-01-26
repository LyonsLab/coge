package CoGe::Factory::PipelineFactory;

use Moose;

use File::Spec::Functions qw(catfile);
use Data::Dumper;

use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Gff;
use CoGe::Builder::Export::Genome;
use CoGe::Builder::Export::Experiment;
use CoGe::Builder::Load::Experiment;
use CoGe::Builder::Load::SRA;
use CoGe::Builder::Load::BatchExperiment;
use CoGe::Builder::Load::Genome;
use CoGe::Builder::Load::Annotation;
use CoGe::Builder::SNP::IdentifySNPs;
use CoGe::Builder::Tools::NucCounter;
use CoGe::Builder::Tools::CoGeBlast;
use CoGe::Builder::Tools::SynMap;
use CoGe::Builder::Tools::SynMap3D;
use CoGe::Builder::Expression::qTeller;
use CoGe::Builder::Methylation::Metaplot;
use CoGe::Builder::PopGen::SummaryStats;
use CoGe::Exception::Generic;

my %typeToClass = (
    'blast'                 => 'CoGe::Builder::Tools::CoGeBlast',
    'export_gff'            => 'CoGe::Builder::Export::Gff',
    'export_fasta'          => 'CoGe::Builder::Export::Fasta',
    'export_genome'         => 'CoGe::Builder::Export::Genome',
    'export_experiment'     => 'CoGe::Builder::Export::Experiment',
    'load_experiment'       => 'CoGe::Builder::Load::Experiment',
    'load_sra'              => 'CoGe::Builder::Load::SRA',
    #'load_batch'            => 'CoGe::Builder::Load::BatchExperiment',
    'load_genome'           => 'CoGe::Builder::Load::Genome',
    'load_annotation'       => 'CoGe::Builder::Load::Annotation',
    'analyze_snps'          => 'CoGe::Builder::SNP::SNPFinder',
    'synmap'                => 'CoGe::Builder::Tools::SynMap',
    'synmap3d'              => 'CoGe::Builder::Tools::SynMap3D',
    'analyze_expression'    => 'CoGe::Builder::Expression::qTeller',
    'analyze_metaplot'      => 'CoGe::Builder::Methylation::Metaplot',
    'analyze_diversity'     => 'CoGe::Builder::PopGen::SummaryStats',
    'analyze_nucleotides'   => 'CoGe::Builder::Tools::NucCounter'
);

sub get {
    my ($self, $request) = @_;

    # Determine pipeline builder
    my $className = $typeToClass{$request->type};
    unless ($className) {
        CoGe::Exception::Generic->throw(message => "Unrecognized job type: " . $request->type);
    }
    my $builder = $className->new(request => $request);

    #
    # Construct the workflow
    #

    # This loosely corresponds to the Extract-Transform-Load paradigm.

    # Pre-build (extract): initialize pipeline and add data retrieval tasks
    my $dr = $builder->pre_build();

    # Build (transform/load): add application-specific pipeline tasks #TODO data passing needs improvement
    $builder->build(
        $dr ? (
                data_files => $dr->data_files,
                data_dir   => $dr->data_dir,
                ncbi_accns => $dr->ncbi_accns
              )
            : ()
    );

    # Post-build: add completion tasks (such as sending notifiation email)
    $builder->post_build();

    return $builder;
}

__PACKAGE__->meta->make_immutable;

1;
