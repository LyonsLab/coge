package CoGe::Factory::PipelineFactory;

use Moose;

use Switch;
use File::Spec::Functions qw(catfile);
use Data::Dumper;

use CoGe::Core::Storage qw(get_workflow_paths);
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
use CoGe::Builder::Tools::CoGeBlast;
use CoGe::Builder::Tools::NucCounter;
use CoGe::Builder::Tools::SynMap;
use CoGe::Builder::Tools::SynMap3D;
use CoGe::Builder::Expression::MeasureExpression;
use CoGe::Builder::Methylation::CreateMetaplot;
use CoGe::Builder::PopGen::MeasureDiversity;

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
    isa => 'CoGe::JEX::Jex',
    is => 'ro',
    required => 1
);

sub get {
    my ($self, $message) = @_;

    my $request = {
        params    => $message->{parameters},
        db        => $self->db,
        user      => $self->user,
        conf      => $self->conf
    };

    # Select pipeline builder
    my $builder;
    switch ($message->{type}) {
        case "blast"              { $builder = CoGe::Builder::Tools::CoGeBlast->new($request) }
        case "export_gff"         { $builder = CoGe::Builder::Export::Gff->new($request) }
        case "export_fasta"       { $builder = CoGe::Builder::Export::Fasta->new($request) }
        case "export_genome"      { $builder = CoGe::Builder::Export::Genome->new($request) }
        case "export_experiment"  { $builder = CoGe::Builder::Export::Experiment->new($request) }
        case "load_experiment"    { $builder = CoGe::Builder::Load::Experiment->new($request) }
        case "load_sra"           { $builder = CoGe::Builder::Load::SRA->new($request) }
        case "load_batch"         { $builder = CoGe::Builder::Load::BatchExperiment->new($request) }
        case "load_genome"        { $builder = CoGe::Builder::Load::Genome->new($request) }
        case "load_annotation"    { $builder = CoGe::Builder::Load::Annotation->new($request) }
        case "analyze_snps"       { $builder = CoGe::Builder::SNP::IdentifySNPs->new($request) }
        case "synmap"             { $builder = CoGe::Builder::Tools::SynMap->new($request) }
        case "synmap3d"           { $builder = CoGe::Builder::Tools::SynMap3D->new($request) }
        case "analyze_expression" { $builder = CoGe::Builder::Expression::MeasureExpression->new($request) }
        case "analyze_metaplot"   { $builder = CoGe::Builder::Methylation::CreateMetaplot->new($request) }
        case "analyze_diversity"  { $builder = CoGe::Builder::PopGen::MeasureDiversity->new($request) }
        case "nuccounter"         { $builder = CoGe::Builder::Tools::NucCounter->new($request) }
        default {
            print STDERR "PipelineFactory::get unknown job type\n";
            return;
        }
    }

    #
    # Construct the workflow
    #

    $builder->pre_build(jex => $self->jex, requester => $message->{requester});
    my $rc = $builder->build();
    unless ($rc) {
        $rc = 'undef' unless defined $rc;
        print STDERR "PipelineFactory::get build failed, rc=$rc\n";
        return;
    }

    # Add completion tasks (such as sending notifiation email)
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
        print $fh Dumper $message, "\n";
        close($fh);
    }
    
    return $builder;
}

1;
