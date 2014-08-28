package CoGe::Pipelines::PipelineFactory;

use Moose;

use CoGe::Pipelines::Misc::Gff;
use CoGe::Accessory::Workflow;

has 'conf' => (
    is => 'ro',
    required => 1
);

sub get {
    my ($self, $message) = @_;

    my ($options, $params) = ($message->{options}, $message->{parameters});
    my $workflow;

    if ($message->{type} eq "gff_export") {
        $workflow = CoGe::Accessory::Workflow->new(name => "generate gff");
        my ($output, %job) = generate_gff($params, $self->conf);
        $workflow->add_job(%job);
    }

    return $workflow;
}

1;
