package CoGe::Builder::GffBuilder;

use Moose;

use CoGe::Core::Storage;
use CoGe::Pipelines::Misc::Gff;
use CoGe::Pipelines::Misc::IPut;

use File::Spec::Functions;
use Data::Dumper;

sub build {
    my $self = shift;

    $self->init_workflow($self->jex);
    return unless $self->workflow->id;

    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name,
                                                        $self->workflow->id);

    my $dest_type = $self->options->{dest_type} or "http";

    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    my ($output, %job) = generate_gff($self->params, $self->conf);
    $self->workflow->add_job(%job);

    if ($dest_type eq "irods") {
        %job = export_to_irods($output, $self->options, $self->user);
        $self->workflow->add_job(%job);
    }

    $self->workflow->add_job(generate_results($output, $dest_type, $result_dir, $self->conf));
}

sub init_workflow {
    my ($self, $jex) = @_;

    $self->workflow($jex->create_workflow(name => "Generate gff", init => 1));
}

with qw(CoGe::Builder::Buildable);

1;
