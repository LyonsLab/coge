package CoGe::Builder::Expression::MeasureExpression;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Core::Storage qw(get_experiment_files);
use CoGe::Builder::Expression::qTeller;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub get_name {
    #my $self = shift;
    return 'Measure Expression'; #TODO add source experiment info
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    unless ($eid) {
        CoGe::Exception::MissingField->throw(message => "Missing experiment_id");
    }
    unless ($self->params->{expression_params}) {
        CoGe::Exception::MissingField->throw(message => "Missing expression_params");
    }
    
    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    unless ($experiment) {
        CoGe::Exception::Generic->throw(message => "Experiment $eid not found");
    }
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

__PACKAGE__->meta->make_immutable;

1;