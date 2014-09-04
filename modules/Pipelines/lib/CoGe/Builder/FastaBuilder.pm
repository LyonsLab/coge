package CoGe::Builder::FastaBuilder;

use Moose;

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Web qw(url_for);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Pipelines::Common::Results;
use CoGe::Pipelines::Misc::IPut;

use File::Spec::Functions;
use Data::Dumper;

sub build {
    my $self = shift;

    $self->init_workflow($self->jex);
    return unless $self->workflow->id;

    my (undef, $result_dir) = get_workflow_paths($self->user->name,
                                                 $self->workflow->id);

    my $dest_type = $self->options->{dest_type};
    $dest_type = "http" unless $dest_type;

    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    my $gid = $self->params->{gid};
    my $genome = get_genome_file($gid);

    if ($dest_type eq "irods") {
        my $base = $self->options->{dest_path};
        $base = irods_get_base_path($self->user->name) unless $base;
        my $dest = catfile($base, "genome_$gid.faa");

        $self->workflow->add_job(export_to_irods($genome, $dest));
        $self->workflow->add_job(generate_results($dest, $dest_type, $result_dir, $self->conf));
    } else {
        $self->workflow->add_job(link_results($genome, $result_dir, $self->conf));
    }
}

sub init_workflow {
    my ($self, $jex) = @_;

    $self->workflow($jex->create_workflow(name => "Get fasta file", init => 1));
}

with qw(CoGe::Builder::Buildable);

1;
