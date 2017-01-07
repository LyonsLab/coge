package CoGe::Builder::PopGen::MeasureDiversity;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Builder::PopGen::SummaryStats qw(build);

sub get_name {
    #my $self = shift;
    return 'Measure Diversity';
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    return unless $eid;
    #return unless $self->params->{diversity_params};
    
    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    return unless $experiment;
    
    #
    # Build workflow steps
    #
    my @tasks;
    
    # Add expression analysis workflow
    my $workflow = CoGe::Builder::PopGen::SummaryStats::build(
        user => $self->user,
        wid => $self->workflow->id,
        experiment => $experiment,
        params => $self->params->{diversity_params}
    );
    push @tasks, @{$workflow->{tasks}};
        
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;