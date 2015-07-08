package CoGe::Builder::Load::Genome;

use Moose;

use Data::Dumper qw(Dumper);
use Switch;
use File::Spec::Functions qw(catfile);

use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path);
use CoGe::Core::Notebook qw(load_notebook);
use CoGe::Builder::CommonTasks;

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

    # Get genome
    my $organism = $self->db->resultset('Organism')->find($organism_id);
    return unless $organism;
    
    # Initialize workflow
    my $info = '"' . $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    $self->workflow($self->jex->create_workflow(name => "Load Genome " . $info, init => 1));
    return unless ($self->workflow && $self->workflow->id);
    
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));

    # Determine file type if not set
    my $file_type = $data->[0]->{file_type}; # type of first data file
    ($file_type) = detect_data_type($file_type, $data->[0]->{path}) unless $file_type;
    
    #
    # Build workflow
    #
    my (@tasks, @input_files, @done_files, @additional_metadata);
    
    # Create tasks to retrieve files #TODO Merge with Load/Experiment.pm
    my $upload_dir = get_upload_path($self->user->name, $load_id);
    foreach my $item (@$data) {
        my $type = lc($item->{type});
        my $task;
        
        # Check if the file already exists which will be the case if called
        # via the LoadGenome page.  
        #TODO Change the LoadGenome page to not actually transfer files in
        # an apache process -- instead let it be done here in the workflow.
        my $filepath = catfile($upload_dir, $item->{path});
        if (-r $filepath) {
            push @input_files, $filepath;
            next;
        }
        
        # Create task based on source type (IRODS, HTTP, FTP)
        if ($type eq 'irods') {
            my $irods_path = $item->{path};
            $irods_path =~ s/^irods//; # strip of leading "irods" from LoadGenome page # FIXME remove this in FileSelect
            $task = create_iget_job(irods_path => $irods_path, local_path => $upload_dir);
            return unless $task;
        }
        elsif ($type eq 'http' or $type eq 'ftp') {
            #TODO
        }
        
        # Add task to workflow
        $self->workflow->add_job($task);
        push @input_files, $task->{outputs}[0];
    }
    
    # Submit workflow to generate genome
    my $job = create_load_genome_job(
        user => $self->user,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        wid => $self->workflow->id,
        gid => $genome->id,
        input_file => $input_files[0],
        metadata => $metadata,
        normalize => $self->params->{normalize} ? $self->params->{normalize_method} : 0
    );
    push @tasks, $job;
    push @done_files, $job->{outputs}->[1];
    
    # Send notification email #TODO merge with Load/Experiment.pm
	if ( $self->params->{email} ) {
	    # Get tiny link
	    my $link;
	    if ($self->requester && $self->requester->{page}) {
	        $link = CoGe::Accessory::Web::get_tiny_link( 
	            url => $self->conf->{SERVER} . $self->requester->{page} . "?wid=" . $self->workflow->id 
	        );
	    }
	    
	    # Build message body
	    my $body = 'Genome "' . $metadata->{name} . '" has finished loading.';
        $body .= "\nLink: $link" if $link;
        $body .= "\n\nNote: you received this email because you submitted a genome on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
	    
		# Create task
		push @tasks, send_email_job(
			from => 'CoGe Support <coge.genome@gmail.com>',
			to => $self->user->email,
			subject => 'CoGe Load Genome done',
			body => $body,
			staging_dir => $staging_dir,
			done_files => \@done_files
		);
	}

    #print STDERR Dumper \@tasks, "\n";
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

with qw(CoGe::Builder::Buildable);

1;
