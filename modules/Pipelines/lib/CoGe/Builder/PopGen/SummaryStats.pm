package CoGe::Builder::PopGen::SummaryStats;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use CoGe::Builder::CommonTasks;

sub get_name {
    #my $self = shift;
    return 'Compute Summary Stats';
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    return unless $eid;
    return unless $self->params->{diversity_params};
    
    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    return unless $experiment;
    my $genome = $experiment->genome;
    # TODO add permissions check here -- or will it happen in Request::Genome?
    
    return 1;
}

1;