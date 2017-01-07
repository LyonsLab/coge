package CoGe::Builder::Expression::MeasureExpression;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Core::Storage qw(get_experiment_files);
use CoGe::Builder::Expression::qTeller qw(build);

sub get_name {
    #my $self = shift;
    return 'Measure Expression';
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    return unless $eid;
    return unless $self->params->{expression_params};
    
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
    
    # Add expression analysis workflow
    my $expression_workflow = CoGe::Builder::Expression::qTeller::build(
        user => $self->user,
        wid => $self->workflow->id,
        genome => $genome,
        input_file => $bam_file,
        metadata => $metadata,
        additional_metadata => $self->params->{additional_metadata},
        params => $self->params->{expression_params}
    );
    push @tasks, @{$expression_workflow->{tasks}};
        
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

1;