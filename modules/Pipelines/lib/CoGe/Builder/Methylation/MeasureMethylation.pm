package CoGe::Builder::Methylation::MeasureMethylation;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use Switch;
use CoGe::Core::Storage qw(get_experiment_files);
use CoGe::Builder::Methylation::Bismark qw(build);
use CoGe::Builder::Methylation::BWAmeth qw(build);

sub get_name {
    #my $self = shift;
    return 'Measure Methylation';
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    return unless $eid;
    return unless $self->params->{methylation_params};
    my $method = $self->params->{methylation_params}->{method};
    
    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    return unless $experiment;
    my $genome = $experiment->genome;
    
    # Copy metadata from input experiment
    my $metadata = { # could almost use experiment->to_hash here except for source_name
        name => $experiment->name,
        version => $experiment->version,
        source => $experiment->source->name,
        restricted => $experiment->restricted
    };
    
    # Get input file
    my $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];

    #
    # Build workflow steps
    #
    my @tasks;
    
    # Add methylation analysis workflow
    my $methyl_params = {
        user => $self->user,
        wid => $self->workflow->id,
        genome => $genome,
        input_file => $bam_file,
        metadata => $metadata,
        read_params => $self->params->{read_params},
        methylation_params => $self->params->{methylation_params}
    };
    
    my $methyl_workflow;
    switch ($method) { # TODO move this into subroutine
        case 'bismark' { $methyl_workflow = CoGe::Builder::Methylation::Bismark::build($methyl_params); }
        case 'bwameth' { $methyl_workflow = CoGe::Builder::Methylation::BWAmeth::build($methyl_params); }
        else           { die "unknown methylation method"; }
    }
    push @tasks, @{$methyl_workflow->{tasks}};
        
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

1;