package CoGe::Builder::Expression::Analyzer;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Core::Storage qw(get_experiment_files);
use CoGe::Builder::Expression::Cufflinks;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my %opts = @_;

    my $bam_file;
    if ($opts{data_files}) {
        ($bam_file) = @{$opts{data_files}};
    }
    else { # for when called from ExperimentView
        my $experiment = $self->request->experiment;
        $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
    }

    unless ($self->params->{metadata}) { # user input experiment's metadat, for when called from ExperimentView
        my $experiment = $self->request->experiment;
        $self->params->{metadata} = { # could almost use experiment->to_hash here except for source_name
            name       => $experiment->name,
            version    => $experiment->version,
            source     => $experiment->source->name,
            restricted => $experiment->restricted
        };
    }

    # Validate inputs (that weren't already checked in Request)
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw( message => "Missing metadata" );
    }

    #
    # Build workflow
    #

    # Add expression analysis tasks -- there is only one method at the moment, Cufflinks
    my $analyzer;
    #my $method = lc($self->params->{expression_params}->{method});
    #if ($method eq 'bismark') {
        $analyzer = CoGe::Builder::Expression::Cufflinks->new($self);
        $analyzer->build(data_files => [$bam_file]);
    #}
    #else {
    #    CoGe::Exception::Generic->throw(message => 'Invalid methylation method');
    #}
    $self->add($analyzer);
}

__PACKAGE__->meta->make_immutable;

1;
