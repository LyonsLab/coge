package CoGe::Builder::Load::Experiment;

use Moose;

use Data::Dumper qw(Dumper);
use Switch;
use File::Spec::Functions qw(catfile);

use CoGe::Core::Notebook qw(load_notebook);
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path);
use CoGe::Builder::CommonTasks;
use CoGe::Builder::Common::Alignment qw(build);
use CoGe::Builder::Expression::qTeller qw(build);
use CoGe::Builder::SNP::CoGeSNPs qw(build);
use CoGe::Builder::SNP::Samtools qw(build);
use CoGe::Builder::SNP::Platypus qw(build);
use CoGe::Builder::SNP::GATK qw(build);

sub build {
    my $self = shift;
    
    # Validate inputs
    my $gid = $self->params->{gid};
    return unless $gid;
    my $data = $self->options->{source_data};
    return unless (defined $data && @$data);
    my $metadata = $self->params->{metadata};
    return unless $metadata;
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get genome
    my $genome = $self->db->resultset('Genome')->find($gid);
    return unless $genome;
    # TODO add permissions check here -- or will it happen in Request::Genome?
    
    # Initialize workflow
    my $info = '"' . $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    $self->workflow($self->jex->create_workflow(name => "Load Experiment " . $info, init => 1));
    return unless $self->workflow->id;
    
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));
    
    # Build workflow steps
    my @tasks;
    my $file_type = $data->[0]->{file_type}; # type of first data file
    
    my @done_files;
    if ( $file_type eq 'fastq' || $file_type eq 'bam' ) {
        my $bam_file;
        my $result_count = 0;
         
        # Align fastq file or take existing bam
        my ($alignment_tasks, $alignment_results);
        if ( $file_type && $file_type eq 'fastq' ) {
            # Add alignment workflow
            ($alignment_tasks, $alignment_results) = CoGe::Builder::Common::Alignment::build(
                user => $self->user,
                wid => $self->workflow->id,
                input_files => $data,
                genome => $genome,
                metadata => $metadata,
                options => $self->options,
                trimming_params => $self->params->{trimming_params},
                alignment_params => $self->params->{alignment_params}
            );
            push @tasks, @$alignment_tasks;
            $bam_file = $alignment_results->{bam_file};
            push @done_files, @{$alignment_results->{done_files}};
            $result_count++;
        }
        elsif ( $file_type && $file_type eq 'bam' ) {
            my $upload_dir = get_upload_path($self->user->name, $self->options->{load_id});
            $bam_file = catfile($upload_dir, $data->[0]->{path});
        
            push @tasks, create_load_bam_job(
                user => $self->user,
                metadata => $metadata,
                staging_dir => $staging_dir,
                result_dir => $result_dir,
                wid => $self->workflow->id,
                gid => $gid,
                bam_file => $bam_file
            );
        }
        else { # error -- should never happen
            die "invalid file type";
        }
        
        # Add expression workflow (if specified)
        my ($expression_tasks, $expression_results);
        if ( $self->params->{expression_params} ) {
            ($expression_tasks, $expression_results) = CoGe::Builder::Expression::qTeller::build(
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_file => $bam_file,
                metadata => $metadata,
                options => $self->options,
                params => $self->params->{expression_params}
            );
            push @tasks, @$expression_tasks;
            push @done_files, @{$expression_results->{done_files}};
            $result_count++;
        }
        
        # Add SNP workflow (if specified)
        my ($snp_tasks, $snp_results);
        if ( $self->params->{snp_params} ) {
            my $method = $self->params->{snp_params}->{method};
            my $snp_params = {
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_file => $bam_file,
                metadata => $metadata,
                options => $self->options,
                params => $self->params->{snp_params},
                skipAnnotations => 1 # annotations for each result experiment are set together in create_notebook_job() later on
            };
            
            switch ($method) { # TODO move this into subroutine
                case 'coge'     { ($snp_tasks, $snp_results) = CoGe::Builder::SNP::CoGeSNPs::build($snp_params); }
                case 'samtools' { ($snp_tasks, $snp_results) = CoGe::Builder::SNP::Samtools::build($snp_params); }
                case 'platypus' { ($snp_tasks, $snp_results) = CoGe::Builder::SNP::Platypus::build($snp_params); }
                case 'gatk'     { ($snp_tasks, $snp_results) = CoGe::Builder::SNP::GATK::build($snp_params); }
                else            { die "unknown SNP method"; }
            }
            push @tasks, @$snp_tasks;
            push @done_files, @{$snp_results->{done_files}};
            $result_count++;
        }
        
        if ($self->options->{notebook}) {
            my @additional_md = (
                qq{https://genomevolution.org/wiki/index.php/Expression_Analysis_Pipeline||note|Generated by CoGe's RNAseq Analysis Pipeline}
            );
            push @additional_md, @{$alignment_results->{metadata}} if ($alignment_results->{metadata});
            push @additional_md, @{$expression_results->{metadata}} if ($expression_results->{metadata});
            push @additional_md, @{$snp_results->{metadata}} if ($snp_results->{metadata});
            
        	push @tasks, add_items_to_notebook_job(
        		notebook_id => $self->options->{notebook_id},
                user => $self->user, 
                wid => $self->workflow->id,
                annotations => \@additional_md,
                staging_dir => $staging_dir,
                done_files => \@done_files        		
        	);
        } 
        else {
	        # Create notebook and additional pipeline metadata depending on results
	        if ( $result_count > 1 ) {
	            my @additional_md = (
	                qq{https://genomevolution.org/wiki/index.php/Expression_Analysis_Pipeline||note|Generated by CoGe's RNAseq Analysis Pipeline}
	            );
	            push @additional_md, @{$alignment_results->{metadata}} if ($alignment_results->{metadata});
	            push @additional_md, @{$expression_results->{metadata}} if ($expression_results->{metadata});
	            push @additional_md, @{$snp_results->{metadata}} if ($snp_results->{metadata});
	            
	            push @tasks, create_notebook_job(
	                user => $self->user,
	                wid => $self->workflow->id,
	                metadata => $metadata, 
	                annotations => \@additional_md,
	                staging_dir => $staging_dir,
	                done_files => \@done_files
	            );
	        }
        }
    }
    # Else, all other file types
    else {
        # Setup full path to input data file
        my $upload_dir = get_upload_path($self->user->name, $self->options->{load_id});
        my $input_file = catfile($upload_dir, $data->[0]->{path});
        
        # Submit workflow to generate experiment
        my $job = create_load_experiment_job(
            user => $self->user,
            staging_dir => $staging_dir,
            result_dir => $result_dir,
            wid => $self->workflow->id,
            gid => $genome->id,
            input_file => $input_file,
            metadata => $metadata,
            normalize => $self->options->{normalize} ? $self->options->{normalize_method} : 0
        );
        push @tasks, $job;
        push @done_files, $job->{outputs}->[1];
    }

	if ( $self->options->{email} ) {
	    # Get tiny link
	    my $link;
	    if ($self->requester && $self->requester->{page}) {
	        $link = CoGe::Accessory::Web::get_tiny_link( 
	            url => $self->conf->{SERVER} . $self->requester->{page} . "?wid=" . $self->workflow->id 
	        );
	    }
	    
	    # Build message body
	    my $body = 'Experiment "' . $metadata->{name} . '" has finished loading.';
        $body .= "\nLink: $link" if $link;
        $body .= "\n\nNote: you received this email because you submitted an experiment on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
	    
		# Create task
		push @tasks, send_email_job(
			from => 'CoGe Support <coge.genome@gmail.com>',
			to => $self->user->email,
			subject => 'CoGe Load Experiment done',
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
