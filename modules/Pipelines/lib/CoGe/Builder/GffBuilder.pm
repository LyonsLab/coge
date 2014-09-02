package CoGe::Builder::GffBuilder;

use Moose;

use CoGe::Pipelines::Misc::Gff;

use File::Spec::Functions;
use Data::Dumper;

sub build {
    my $self = shift;

    $self->init_workflow($self->jex);
    return unless $self->workflow->id;

    my $result_dir = catdir($self->conf->{SECTEMPDIR}, "results", $self->user->name, $self->workflow->id);
    my $dest_type = $self->options->{dest_type};
    $dest_type = "http" unless $dest_type;

    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    my ($output, %job) = generate_gff($self->params, $self->conf);
    $self->workflow->add_job(%job);
    $self->workflow->add_job(generate_results($output, $dest_type, $result_dir, $self->conf));
}

sub init_workflow {
    my ($self, $jex) = @_;

    $self->workflow($jex->create_workflow(name => "Generate gff", init => 1));
}

with qw(CoGe::Builder::Buildable);

1;
