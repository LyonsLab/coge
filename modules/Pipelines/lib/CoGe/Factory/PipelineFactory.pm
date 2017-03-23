package CoGe::Factory::PipelineFactory;

use Moose;

use File::Spec::Functions qw(catfile);

use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Gff;
use CoGe::Builder::Export::Genome;
use CoGe::Builder::Export::Experiment;
use CoGe::Builder::Load::Experiment;
use CoGe::Builder::Load::SRA;
use CoGe::Builder::Load::BatchExperiment;
use CoGe::Builder::Load::Genome;
use CoGe::Builder::Load::Annotation;
use CoGe::Builder::Tools::PercentGCAT;
use CoGe::Builder::Tools::CoGeBlast;
use CoGe::Builder::Tools::SynMap;
use CoGe::Builder::Tools::SynMap3D;
use CoGe::Builder::Tools::Pseudoassembly;
use CoGe::Builder::SNP::Analyzer;
use CoGe::Builder::SNP::Merge;
use CoGe::Builder::Expression::Analyzer;
use CoGe::Builder::Methylation::Analyzer;
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
    'synmap'                => 'CoGe::Builder::Tools::SynMap',
    'synmap3d'              => 'CoGe::Builder::Tools::SynMap3D',
    'pseudoassembly'        => 'CoGe::Builder::Tools::Pseudoassembly',
    'analyze_snps'          => 'CoGe::Builder::SNP::Analyzer',
    'analyze_expression'    => 'CoGe::Builder::Expression::Analyzer',
    'analyze_methylation'   => 'CoGe::Builder::Methylation::Analyzer',
    'analyze_metaplot'      => 'CoGe::Builder::Methylation::Metaplot',
    'analyze_diversity'     => 'CoGe::Builder::PopGen::SummaryStats',
    'merge_snps'            => 'CoGe::Builder::SNP::Merge',
    'percent_gc_at'         => 'CoGe::Builder::Tools::PercentGCAT'
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
