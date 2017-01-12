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
#    my $self = shift;
#
#    # Validate inputs
#    my $eid = $self->params->{eid} || $self->params->{experiment_id};
#    unless ($eid) {
#        CoGe::Exception::MissingField->throw(message => "Missing experiment_id");
#    }
#    unless ($self->params->{expression_params}) {
#        CoGe::Exception::MissingField->throw(message => "Missing expression_params");
#    }
#
#    # Get experiment
#    my $experiment = $self->request->experiment;
#    my $genome = $self-request->genome;
#
#    # Copy metadata from input experiment
#    my $metadata = { # could almost use experiment->to_hash here except for source_name
#        name => $experiment->name,
#        version => $experiment->version,
#        source => $experiment->source->name,
#        restricted => $experiment->restricted
#    };
#
#    # Get input file
#    my $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
#
#    # Add expression analysis tasks
#    my $expression_workflow = CoGe::Builder::Expression::qTeller::new->(
#    );
#    $expr->build($bam_file);
}

__PACKAGE__->meta->make_immutable;

1;