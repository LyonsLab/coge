package CoGe::Builder::GenomeBuilder;

use Moose;

use CoGe::Builder::FastaBuilder;
use CoGe::Builder::GffBuilder;
use CoGe::Core::Storage;
use CoGe::Accessory::Utils;
use CoGe::Pipelines::Misc::Gff;
use CoGe::Pipelines::Misc::IPut;
use CoGe::Pipelines::Common::Results;

use File::Spec::Functions qw(catfile catdir);
use Data::Dumper;

sub build {
    my $self = shift;

    my $dest_type = $self->options->{dest_type};
       $dest_type = "http" unless $dest_type;

    # Initialize workflow
    $self->workflow( $self->jex->create_workflow(name => "Export genome", init => 1) );
    return unless $self->workflow->id;
    
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));
    
    #
    # FASTA
    #
    
    # Get genome data file path
    return unless (defined $self->params && $self->params->{gid});
    my $gid = $self->params->{gid};
    my $genome = $self->db->resultset("Genome")->find($gid);
    my $genome_file = get_genome_file($gid);

    # Determine name of exported file
    my $genome_name = sanitize_name($genome->organism->name);#escape($genome->organism->name);
       $genome_name = 'genome_'.$gid unless $genome_name;
    my $output_file = $genome_name.'.faa';

    # Setup tasks to export/download the file
#    if ($dest_type eq "irods") { # irods export
#        my $irods_base = $self->options->{dest_path};
#        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
#
#        my $irods_dest = catfile($irods_base, $output_file);
#        my $irods_done = catfile($staging_dir, "irods.done");
#
#        $self->workflow->add_job( export_to_irods($genome_file, $irods_dest, $self->options->{overwrite}, $irods_done) );
#        $self->workflow->add_job( generate_results($irods_dest, $dest_type, $result_dir, $self->conf, $irods_done) );
#    }
#    else { # http download
#        $self->workflow->add_job( link_results($genome_file, $output_file, $result_dir, $self->conf) );
#    }
    
    #
    # GFF
    #
    
    $self->params->{basename} = sanitize_name($genome_name);
    my ($output_gff_file, %job) = generate_gff($self->params, $self->conf);
    $self->workflow->add_job(%job);

#    if ($dest_type eq "irods") { # irods export
#        my $irods_base = $self->options->{dest_path};
#        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
#        my $irods_dest = catfile($irods_base, basename($output));
#        my $irods_done = catfile($staging_dir, "irods.done");
#
#        $self->workflow->add_job( export_to_irods($output, $irods_dest, $self->options->{overwrite}, $irods_done) );
#        $self->workflow->add_job( generate_results($irods_dest, $dest_type, $result_dir, $self->conf, $irods_done) );
#    } 
#    else { # http download
#        $self->workflow->add_job( link_results($output, $output, $result_dir, $self->conf) );
#    }
    
    #
    # Create tarball
    #
    
    my $filename = 'genome_'.$genome_name.'.tar.gz';
    my $output_dir = $self->get_download_path($gid);
    my $output_archive_file = catfile($output_dir, $filename);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    %job = export_genome($gid, $output_archive_file, $output_gff_file, $self->conf);
    $self->workflow->add_job(%job);

    if ($dest_type eq "irods") {
        my $base = $self->options->{dest_path};
        $base = irods_get_base_path($self->user->name) unless $base;
        my $dest = catfile($base, $filename);
        my $irods_done = catfile($staging_dir, "irods.done");

        $self->workflow->add_job( export_to_irods($output_archive_file, $dest, $self->options->{overwrite}, $irods_done) );
        $self->workflow->add_job( generate_results($dest, $dest_type, $result_dir, $self->conf, $irods_done) );
    } 
    else {
        $self->workflow->add_job( link_results($output_archive_file, $output_archive_file, $result_dir, $self->conf) );
    }
    
    return 1;
}

sub get_download_path { #TODO merge with similar subroutine in ExperimentBuilder.pm and move into Storage.pm
    my $self = shift;
    my $gid = shift;
    my $unique_path = get_unique_id();
    my @paths = ($self->conf->{SECTEMPDIR}, 'downloads', 'genomes', $gid, $unique_path);
    return File::Spec->catdir(@paths);
}

sub export_genome {
    my ($gid, $output_file, $gff_file, $conf) = @_;

    return (
        cmd => catdir($conf->{SCRIPTDIR}, "export_experiment_or_genome.pl"),
        description => "Generating genome files",
        args => [
            ["-id", $gid, 0],
            ["-type", 'genome', 0],
            ["-output", $output_file, 1],
            ["-conf", $conf->{_CONFIG_PATH}, 0],
            ["-dir", ".", ""],
            ["-files", $gff_file, 0]
        ],
        inputs => [],
        outputs => [$output_file]
    );
}

with qw(CoGe::Builder::Buildable);

1;
