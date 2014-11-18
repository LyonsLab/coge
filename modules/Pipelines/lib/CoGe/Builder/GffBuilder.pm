package CoGe::Builder::GffBuilder;

use Moose;

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Core::Storage;
use CoGe::Pipelines::Common::Results;
use CoGe::Pipelines::Misc::Gff;
use CoGe::Pipelines::Misc::IPut;

use File::Basename qw(basename);
use File::Spec::Functions;
use Data::Dumper;

sub build {
    my $self = shift;

    $self->init_workflow($self->jex);
    return unless $self->workflow->id;

    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);

    my $dest_type = $self->options->{dest_type};
    $dest_type = "http" unless $dest_type;

    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    my $genome = $self->db->resultset("Genome")->find($self->params->{gid});
    $self->params->{basename} = sanitize_name($genome->organism->name);

    my ($output, %job) = generate_gff($self->params, $self->conf);
    $self->workflow->add_job(%job);

    if ($dest_type eq "irods") { # irods export
        my $irods_base = $self->options->{dest_path};
        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
        my $irods_dest = catfile($irods_base, basename($output));
        my $irods_done = catfile($staging_dir, "irods.done");

        $self->workflow->add_job( export_to_irods($output, $irods_dest, $self->options->{overwrite}, $irods_done) );
        $self->workflow->add_job( generate_results($irods_dest, $dest_type, $result_dir, $self->conf, $irods_done) );
    } 
    else { # http download
        $self->workflow->add_job( link_results($output, $output, $result_dir, $self->conf) );
    }
}

sub init_workflow {
    my ($self, $jex) = @_;

    $self->workflow($jex->create_workflow(name => "Generate gff", init => 1));
}

with qw(CoGe::Builder::Buildable);

1;
