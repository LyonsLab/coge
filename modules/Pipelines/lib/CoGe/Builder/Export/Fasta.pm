package CoGe::Builder::Export::Fasta;

use Moose;
with qw(CoGe::Builder::Buildable);

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Accessory::TDS;
use CoGe::Core::Storage qw(get_genome_file);
use CoGe::Core::Genome qw(get_irods_metadata);
use CoGe::Builder::CommonTasks qw(generate_gff export_to_irods generate_results link_results send_email_job create_irods_imeta_job);

use File::Spec::Functions qw(catfile);
use Data::Dumper;
#use URI::Escape::JavaScript qw(escape);

sub get_name {
    return "Export FASTA";
}

sub build {
    my $self = shift;

    # Get genome data file path
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
    my @done_files;
    my $dest_type = $self->params->{dest_type} // 'http';
    if ($dest_type eq "irods") { # irods export
        # Set IRODS destination path
        my $irods_base = $self->params->{dest_path};
        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
        my $irods_dest = catfile($irods_base, $output_file);

        # Export file task (IRODS iput)
        my $irods_done = catfile($self->staging_dir, "irods.done");
        $self->workflow->add_job( export_to_irods($genome_file, $irods_dest, $self->params->{overwrite}, $irods_done) );
        
        # Set file metadata task (IRODS imeta)
        my $md = get_irods_metadata($genome);
        my $md_file = catfile(get_genome_cache_path($genome->id), 'irods_metadata.json');
        TDS::write($md_file, $md);
        $self->workflow->add_job(
            create_irods_imeta_job(
                dest_file_path => $irods_dest,
                metadata => $md_file
            )
        );
        
        #TODO remove legacy generate_results
        my $results_task = generate_results($irods_dest, $dest_type, $self->result_dir, $self->conf, $irods_done);
        $self->workflow->add_job($results_task);
        push @done_files, $results_task->{outputs}->[0];
    }
    else { # http download
        my $results_task = link_results($genome_file, $output_file, $self->result_dir, $self->conf);
        $self->workflow->add_job($results_task);
        push @done_files, $results_task->{outputs}->[0];
    }
    
    # Send notification email #TODO move into shared module
    if ( $self->params->{email} ) {
        # Build message body
        my $body = 'FASTA export for genome "' . $genome_name . '" (id' . $gid . ') has finished.';
        $body .= "\nLink: " . $self->site_url if $self->site_url;
        $body .= "\n\nNote: you received this email because you submitted a job on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
        
        # Create task
        $self->workflow->add_job(
            send_email_job(
                to => $self->user->email,
                subject => 'CoGe FASTA export done',
                body => $body,
                staging_dir => $self->staging_dir,
                done_files => \@done_files
            )
        );
    }
    
    return 1;
}

1;
