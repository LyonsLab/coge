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
use CoGe::Builder::Expression::MeasureExpression;
use CoGe::Builder::Methylation::CreateMetaplot;
use CoGe::Builder::PopGen::MeasureDiversity;

my %typeToClass = (
    'blast'                 => 'CoGe::Builder::Tools::CoGeBlast',
    'export_gff'            => 'CoGe::Builder::Export::Gff',
    'export_fasta'          => 'CoGe::Builder::Export::Fasta',
    'export_genome'         => 'CoGe::Builder::Export::Genome',
    'export_experiment'     => 'CoGe::Builder::Export::Experiment',
    'load_experiment'       => 'CoGe::Builder::Load::Experiment',
    'load_sra'              => 'CoGe::Builder::Load::SRA',
    'load_batch'            => 'CoGe::Builder::Load::BatchExperiment',
    'load_genome'           => 'CoGe::Builder::Load::Genome',
    'load_annotation'       => 'CoGe::Builder::Load::Annotation',
    'analyze_snps'          => 'CoGe::Builder::SNP::IdentifySNPs',
    'synmap'                => 'CoGe::Builder::Tools::SynMap',
    'synmap3d'              => 'CoGe::Builder::Tools::SynMap3D',
    'analyze_expression'    => 'CoGe::Builder::Expression::MeasureExpression',
    'analyze_metaplot'      => 'CoGe::Builder::Methylation::CreateMetaplot',
    'analyze_diversity'     => 'CoGe::Builder::PopGen::MeasureDiversity',
    'nuccounter'            => 'CoGe::Builder::Tools::NucCounter'
);

sub get {
    my ($self, $request) = @_;

    # Determine pipeline builder
    my $className = $typeToClass{$request->type};
    unless ($className) {
        warn "PipelineFactory::get unknown job type";
        return;
    }

    my $builder = $className->new(request => $request);

    #
    # Construct the workflow
    #

    # Pre-build: initialize pipeline
    $builder->pre_build();

    # Build: add application-specific pipeline tasks
    my $rc = $builder->build();
    unless ($rc) {
        $rc = 'undef' unless defined $rc;
        warn "PipelineFactory::get build failed, rc=$rc";
        return;
    }

    # Post-build: add completion tasks (such as sending notifiation email)
    $builder->post_build();

    # Dump info to file for debugging
    if ($builder->result_dir) {
        my $cmd = 'chmod g+rw ' . $builder->result_dir;
        `$cmd`;

        # Dump raw workflow
        open(my $fh, '>', catfile($builder->result_dir, 'workflow.log'));
        print $fh Dumper $builder->workflow, "\n";
        close($fh);

        # Dump params
        open($fh, '>', catfile($builder->result_dir, 'params.log'));
        print $fh Dumper $request->payload, "\n";
        close($fh);
    }
    
    return $builder;
}

__PACKAGE__->meta->make_immutable;

1;
