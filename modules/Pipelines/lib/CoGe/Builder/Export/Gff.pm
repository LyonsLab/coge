package CoGe::Builder::Export::Gff;

use Moose;
with qw(CoGe::Builder::Buildable);

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Builder::CommonTasks qw(generate_gff export_to_irods generate_results link_results send_email_job);

use File::Basename qw(basename);
use File::Spec::Functions;
use Data::Dumper;

sub get_name {
    return "Generate/export gff";
}

sub build {
    my $self = shift;

    # Verify required parameters and set defaults
    my $dest_type = $self->params->{dest_type};
    $dest_type = "http" unless $dest_type;
    
    my $gid = $self->params->{gid} || $self->params->{genome_id};
    $self->params->{gid} = $gid; # required by generate_gff() call below
    return unless $gid;

    # Get genome
    my $genome = $self->db->resultset("Genome")->find($gid);
    my $genome_name = $self->params->{basename} = sanitize_name($genome->organism->name);
    $genome_name = 'genome_'.$gid unless $genome_name;

    # Generate GFF file
    my @done_files;
    my ($output, $task) = generate_gff(%{$self->params});
    $self->workflow->add_job($task);

    if ($dest_type eq "irods") { # irods export
        my $irods_base = $self->params->{dest_path};
        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
        my $irods_dest = catfile($irods_base, basename($output));
        my $irods_done = catfile($self->staging_dir, "irods.done");

        $self->workflow->add_job( export_to_irods($output, $irods_dest, $self->params->{overwrite}, $irods_done) );
        my $results_task = generate_results($irods_dest, $dest_type, $self->result_dir, $self->conf, $irods_done);
        $self->workflow->add_job($results_task);
        push @done_files, $results_task->{outputs}->[0];
    } 
    else { # http download
        my $results_task = link_results($output, $output, $self->result_dir, $self->conf);
        $self->workflow->add_job($results_task);
        push @done_files, $results_task->{outputs}->[0];
    }
    
    # Send notification email #TODO move into shared module
    if ( $self->params->{email} ) {
        # Build message body
        my $body = 'GFF export for genome "' . $genome_name . '" (id' . $gid . ') has finished.';
        $body .= "\nLink: " . $self->site_url if $self->site_url;
        $body .= "\n\nNote: you received this email because you submitted a job on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
        
        # Create task
        $self->workflow->add_job(
            send_email_job(
                to => $self->user->email,
                subject => 'CoGe GFF export done',
                body => $body,
                staging_dir => $self->staging_dir,
                done_files => \@done_files
            )
        );
    }
    
    return 1;
}

1;
