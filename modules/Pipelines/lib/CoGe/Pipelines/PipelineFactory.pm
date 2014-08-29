package CoGe::Pipelines::PipelineFactory;

use Moose;

use File::Spec::Functions;

use CoGe::Pipelines::Misc::Gff;
use CoGe::Accessory::Workflow;

has 'conf' => (
    is => 'ro',
    required => 1
);


has 'user' => (
    is  => 'ro',
    required => 1
);

has 'jex' => (
    is => 'ro',
    required => 1
);

sub get {
    my ($self, $message) = @_;

    my ($options, $params) = ($message->{options}, $message->{parameters});
    my $workflow;

    if ($message->{type} eq "gff_export") {
        $workflow = $self->jex->create_workflow(name => "generate gff", init => 1);

        my $result_dir = catdir($self->conf->{SECTEMPDIR}, "results", $self->user->name, $workflow->id);
        my $dest_type = $options->{dest_type};
        $dest_type = "http" unless $dest_type;

        my ($output, %job) = generate_gff($params, $self->conf);
        $workflow->logfile(catfile($result_dir, "debug.log"));
        $workflow->add_job(%job);
        $workflow->add_job(generate_results($output, $dest_type, $result_dir, $self->conf));
    }

    return $workflow;
}

1;
