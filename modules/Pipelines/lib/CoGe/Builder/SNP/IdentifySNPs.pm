package CoGe::Builder::SNP::IdentifySNPs;

use Moose;

use Data::Dumper qw(Dumper);
use Switch;
use File::Spec::Functions qw(catfile);

use CoGe::Core::Storage qw(get_workflow_paths get_upload_path get_experiment_files);
use CoGe::Builder::CommonTasks;
use CoGe::Builder::SNP::CoGeSNPs qw(build);
use CoGe::Builder::SNP::Samtools qw(build);
use CoGe::Builder::SNP::Platypus qw(build);
use CoGe::Builder::SNP::GATK qw(build);

sub build {
    my $self = shift;
    
    # Validate inputs
    my $eid = $self->params->{eid};
    return unless $eid;
    return unless $self->params->{snp_params};
    my $method = $self->params->{snp_params}->{method};
    
    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    return unless $experiment;
    my $genome = $experiment->genome;
    # TODO add permissions check here -- or will it happen in Request::Genome?
    
    # Copy metadata from input experiment
    my $metadata = { # could almost use experiment->to_hash here except for source_name
        name => $experiment->name,
        version => $experiment->version,
        source => $experiment->source->name,
        restricted => $experiment->restricted
    };
    
    # Get input file
    my $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
    
    # Initialize workflow
    my $workflow_name = 'Identify SNPs' . ( $method ? ' using '.$method.' method' : '' );
    $self->workflow($self->jex->create_workflow(name => $workflow_name, init => 1));
    return unless ($self->workflow && $self->workflow->id);
    
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));
    
    # Build workflow steps
    my @tasks;
    
    # Add SNP workflow
    my $snp_params = {
        user => $self->user,
        wid => $self->workflow->id,
        genome => $genome,
        input_file => $bam_file,
        metadata => $metadata,
        params => $self->params->{snp_params}
    };
    
    my $snp_workflow;
    switch ($method) { # TODO move this into subroutine
        case 'coge'     { $snp_workflow = CoGe::Builder::SNP::CoGeSNPs::build($snp_params); }
        case 'samtools' { $snp_workflow = CoGe::Builder::SNP::Samtools::build($snp_params); }
        case 'platypus' { $snp_workflow = CoGe::Builder::SNP::Platypus::build($snp_params); }
        case 'gatk'     { $snp_workflow = CoGe::Builder::SNP::GATK::build($snp_params); }
        else            { die "unknown SNP method"; }
    }
    push @tasks, @{$snp_workflow->{tasks}};
        
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

with qw(CoGe::Builder::Buildable);

1;