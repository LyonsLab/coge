package CoGe::Builder::Export::Fasta;

use Moose;

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Builder::CommonTasks;

use File::Spec::Functions;
use Data::Dumper;
#use URI::Escape::JavaScript qw(escape);

sub build {
    my $self = shift;

    # Initialize workflow
    $self->workflow( $self->jex->create_workflow(name => "Export FASTA", init => 1) );
    return unless ($self->workflow && $self->workflow->id);

    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);

    my $dest_type = $self->params->{dest_type};
       $dest_type = "http" unless $dest_type;

    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    # Get genome data file path
    return unless (defined $self->params && ($self->params->{gid} || $self->params->{genome_id}));
    my $gid = $self->params->{gid} || $self->params->{genome_id};
    return unless $gid;
    my $genome = $self->db->resultset("Genome")->find($gid);
    return unless $genome;
    my $genome_file = get_genome_file($gid);

    # Determine name of exported file
    my $genome_name = sanitize_name($genome->organism->name);#escape($genome->organism->name);
       $genome_name = 'genome_'.$gid unless $genome_name;
    my $output_file = $genome_name.'.faa';

    # Setup tasks to export/download the file
    if ($dest_type eq "irods") { # irods export
        my $irods_base = $self->params->{dest_path};
        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;

        my $irods_dest = catfile($irods_base, $output_file);
        my $irods_done = catfile($staging_dir, "irods.done");

        $self->workflow->add_job( export_to_irods($genome_file, $irods_dest, $self->params->{overwrite}, $irods_done) );
        $self->workflow->add_job( generate_results($irods_dest, $dest_type, $result_dir, $self->conf, $irods_done) );
    }
    else { # http download
        $self->workflow->add_job( link_results($genome_file, $output_file, $result_dir, $self->conf) );
    }
    
    return 1;
}

with qw(CoGe::Builder::Buildable);

1;
