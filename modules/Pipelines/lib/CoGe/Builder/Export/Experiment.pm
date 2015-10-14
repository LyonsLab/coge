package CoGe::Builder::Export::Experiment;

use Moose;
with qw(CoGe::Builder::Buildable);

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Core::Storage qw(get_download_path);
use CoGe::Builder::CommonTasks qw(export_experiment_job export_to_irods generate_results link_results send_email_job);
use File::Spec::Functions qw(catdir catfile);
use Data::Dumper;

sub get_name {
    return "Export experiment";
}

sub build {
    my $self = shift;

    # Verify required parameters and set defaults
    my $dest_type = $self->params->{dest_type};
    $dest_type = "http" unless $dest_type;
    
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    return unless $eid;

    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    my $exp_name = sanitize_name($experiment->name);
       $exp_name = $eid unless $exp_name;

    my $filename = "experiment_$exp_name.tar.gz";
    my $cache_dir = get_download_path('experiment', $eid);
    my $cache_file = catfile($cache_dir, $filename);

    # Export experiment
    my @done_files;
    $self->workflow->add_job( export_experiment_job(eid => $eid, output => $cache_file) );

    if ($dest_type eq "irods") {
        my $base = $self->params->{dest_path};
        $base = irods_get_base_path($self->user->name) unless $base;
        my $dest = catfile($base, $filename);
        my $irods_done = catfile($self->staging_dir, "irods.done");

        $self->workflow->add_job( export_to_irods($cache_file, $dest, $self->params->{overwrite}, $irods_done) );
        my $results_task = generate_results($dest, $dest_type, $self->result_dir, $self->conf, $irods_done);
        $self->workflow->add_job($results_task);
        push @done_files, $results_task->{outputs}->[0];
    } 
    else {
        my $results_task = link_results($cache_file, $cache_file, $self->result_dir, $self->conf);
        $self->workflow->add_job($results_task);
        push @done_files, $results_task->{outputs}->[0];
    }
    
    # Send notification email #TODO move into shared module
    if ( $self->params->{email} ) {
        # Build message body
        my $body = 'Export of experiment "' . $exp_name . '" has finished.';
        $body .= "\nLink: " . $self->site_url if $self->site_url;
        $body .= "\n\nNote: you received this email because you submitted a job on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
        
        # Create task
        $self->workflow->add_job(
            send_email_job(
                to => $self->user->email,
                subject => 'CoGe experiment export done',
                body => $body,
                staging_dir => $self->staging_dir,
                done_files => \@done_files
            )
        );
    }
    
    return 1;
}

1;
