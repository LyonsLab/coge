package CoGe::Builder::Export::Gff;

use Moose;

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Core::Storage;
use CoGe::Builder::CommonTasks;

use File::Basename qw(basename);
use File::Spec::Functions;
use Data::Dumper;

sub build {
    my $self = shift;

    # Initialize workflow
    $self->workflow($self->jex->create_workflow(name => "Generate/export gff", init => 1));
    return unless ($self->workflow && $self->workflow->id);

    # Verify required parameters and set defaults
    my $dest_type = $self->params->{dest_type};
    $dest_type = "http" unless $dest_type;

    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    # Get genome
    return unless (defined $self->params && $self->params->{gid});
    my $genome = $self->db->resultset("Genome")->find($self->params->{gid});
    $self->params->{basename} = sanitize_name($genome->organism->name);

    my ($output, $task) = generate_gff(%{$self->params});
    $self->workflow->add_job($task);

    if ($dest_type eq "irods") { # irods export
        my $irods_base = $self->params->{dest_path};
        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
        my $irods_dest = catfile($irods_base, basename($output));
        my $irods_done = catfile($staging_dir, "irods.done");

        $self->workflow->add_job( export_to_irods($output, $irods_dest, $self->params->{overwrite}, $irods_done) );
        $self->workflow->add_job( generate_results($irods_dest, $dest_type, $result_dir, $self->conf, $irods_done) );
    } 
    else { # http download
        $self->workflow->add_job( link_results($output, $output, $result_dir, $self->conf) );
    }
    
    return 1;
}

with qw(CoGe::Builder::Buildable);

1;
