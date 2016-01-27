package CoGe::Builder::Load::Genome;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile);

use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_upload_path);
use CoGe::Builder::CommonTasks;

sub get_name {
    my $self = shift;
    my $metadata = $self->params->{metadata};
    my $info;
    $info .= $metadata->{organism} if $metadata->{organism};
    $info .= " (" . $metadata->{name} . ")"  if $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    return "Load Genome \"$info\"";
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $organism_id = $self->params->{organism_id};
    return unless $organism_id;
    my $data = $self->params->{source_data};
    return unless (defined $data && @$data);
    my $metadata = $self->params->{metadata};
    return unless $metadata;
    my $load_id = $self->params->{load_id} || get_unique_id();
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get organism
    my $organism = $self->db->resultset('Organism')->find($organism_id);
    return unless $organism;
    
    #
    # Build workflow
    #
    my (@tasks, @input_files, @done_files, @ncbi_accns);
    
    # Create tasks to retrieve files
    my $upload_dir = get_upload_path($self->user->name, $load_id);
    my $data_workflow = create_data_retrieval_workflow(upload_dir => $upload_dir, data => $data);
    push @tasks, @{$data_workflow->{tasks}} if ($data_workflow->{tasks});
    push @input_files, @{$data_workflow->{outputs}} if ($data_workflow->{outputs});
    push @ncbi_accns, @{$data_workflow->{ncbi}} if ($data_workflow->{ncbi});
    
    # Submit workflow to add genome
    if (@ncbi_accns) { # NCBI-based load
        my $task = create_load_genome_from_NCBI_job(
            user => $self->user,
            staging_dir => $self->staging_dir,
            wid => $self->workflow->id,
            ncbi_accns => \@ncbi_accns,
            metadata => $metadata,
        );
        push @tasks, $task;
        push @done_files, $task->{outputs}->[1];
    }
    else { # File-based load
        my $task = create_load_genome_job(
            user => $self->user,
            staging_dir => $self->staging_dir,
            wid => $self->workflow->id,
            organism_id => $organism->id,
            input_files => \@input_files,
            metadata => $metadata,
        );
        push @tasks, $task;
        push @done_files, $task->{outputs}->[1];
    }
    
    # Send notification email #TODO move into shared module
	if ( $self->params->{email} ) {
	    # Build message body
	    my $body = 'Genome "' . $metadata->{name} . '" has finished loading.';
        $body .= "\nLink: " . $self->site_url if $self->site_url;
        $body .= "\n\nNote: you received this email because you submitted a job on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
	    
		# Create task
		push @tasks, send_email_job(
			to => $self->user->email,
			subject => 'CoGe Load Genome done',
			body => $body,
			staging_dir => $self->staging_dir,
			done_files => \@done_files
		);
	}

    #print STDERR Dumper \@tasks, "\n";
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

1;
