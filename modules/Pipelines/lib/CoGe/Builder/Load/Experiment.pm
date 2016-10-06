package CoGe::Builder::Load::Experiment;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use Switch;
use File::Spec::Functions qw(catfile);

use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path);
use CoGe::Core::Experiment qw(detect_data_type);
use CoGe::Core::Metadata qw(to_annotations);
use CoGe::Builder::CommonTasks;
use CoGe::Builder::Common::Alignment qw(build);
use CoGe::Builder::Expression::qTeller qw(build);
use CoGe::Builder::PopGen::SummaryStats qw(build);
use CoGe::Builder::SNP::CoGeSNPs qw(build);
use CoGe::Builder::SNP::Samtools qw(build);
use CoGe::Builder::SNP::Platypus qw(build);
use CoGe::Builder::SNP::GATK qw(build);
use CoGe::Builder::Methylation::Bismark qw(build);
use CoGe::Builder::Methylation::BWAmeth qw(build);
use CoGe::Builder::Methylation::Metaplot qw(build);
use CoGe::Builder::Protein::ChIPseq qw(build);

sub get_name {
    my $self = shift;
    my $metadata = $self->params->{metadata};
    my $info = '"' . $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    return "Load Experiment " . $info;
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $gid = $self->params->{genome_id};
    return unless $gid;
    my $data = $self->params->{source_data};
    return unless (defined $data && @$data);
    my $metadata = $self->params->{metadata};
    return unless $metadata;
    my $additional_metadata = $self->params->{additional_metadata}; # optional
    my $load_id = $self->params->{load_id} || get_unique_id();
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;
    
    # Get genome
    my $genome = $self->db->resultset('Genome')->find($gid);
    return unless $genome;
    
    # Determine file type if not set
    my $file_type = $data->[0]->{file_type}; # type of first data file
    ($file_type) = detect_data_type($file_type, $data->[0]->{path}) unless $file_type;
    
    #
    # Build workflow
    #
    my (@tasks, @input_files, @done_files);
    
    # Create tasks to retrieve files
    my $upload_dir = get_upload_path($self->user->name, $load_id);
    my $data_workflow = create_data_retrieval_workflow(upload_dir => $upload_dir, data => $data);
    push @tasks, @{$data_workflow->{tasks}} if ($data_workflow->{tasks});
    push @input_files, @{$data_workflow->{outputs}} if ($data_workflow->{outputs});
    
    # Build analytical tasks based on file type
    if ( $file_type eq 'fastq' || $file_type eq 'bam' || $file_type eq 'sra' ) {
        my ($bam_file, $bam_files); #TODO reconcile these
        my @raw_bam_files; # mdb added 2/29/16 for Bismark, COGE-706
         
        # Align fastq file or take existing bam
        if ( $file_type && ( $file_type eq 'fastq' || $file_type eq 'sra' ) ) {
            # Add alignment workflow
            my $alignment_workflow = CoGe::Builder::Common::Alignment::build(
                user => $self->user,
                wid => $self->workflow->id,
                input_files => \@input_files,
                genome => $genome,
                metadata => $metadata,
                additional_metadata => $additional_metadata,
                load_id => $load_id,
                params => $self->params
            );
            return if ($alignment_workflow->{error}); #TODO need to propagate this error up to client
            
            push @tasks, @{$alignment_workflow->{tasks}};
            $bam_files = $alignment_workflow->{bam_files};
            $bam_file = $bam_files->[0];
            @raw_bam_files = @{$alignment_workflow->{raw_bam_files}};
            push @done_files, @{$alignment_workflow->{done_files}};
        }
        elsif ( $file_type && $file_type eq 'bam' ) {
            $bam_file = $input_files[0];
            $bam_files = \@input_files;
            @raw_bam_files = @$bam_files;
            
            my $annotations = CoGe::Core::Metadata::to_annotations($additional_metadata);
            
            my $bam_task = create_load_bam_job(
                user => $self->user,
                metadata => $metadata,
                annotations => $annotations,
                staging_dir => $self->staging_dir,
                wid => $self->workflow->id,
                gid => $gid,
                bam_file => $bam_file
            );
            push @tasks, $bam_task;
            push @done_files, $bam_task->{outputs}->[1];
        }
        else { # error -- should never happen
            die "invalid file type";
        }
        
        # Add expression workflow (if specified)
        if ( $self->params->{expression_params} ) {
            my $expression_workflow = CoGe::Builder::Expression::qTeller::build(
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_file => $bam_file,
                metadata => $metadata,
                additional_metadata => $additional_metadata,
                params => $self->params->{expression_params}
            );
            push @tasks, @{$expression_workflow->{tasks}};
            push @done_files, @{$expression_workflow->{done_files}};
        }
        
        # Add SNP workflow (if specified)
        if ( $self->params->{snp_params} ) {
            my $method = $self->params->{snp_params}->{method};
            my $snp_params = {
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_file => $bam_file,
                sorted => ($file_type ne 'bam'),
                metadata => $metadata,
                additional_metadata => $additional_metadata,
                params => $self->params->{snp_params},
            };
            
            my $snp_workflow;
            switch ($method) { #FIXME pass into IdentifySNPs instead
                case 'coge'     { $snp_workflow = CoGe::Builder::SNP::CoGeSNPs::build($snp_params); }
                case 'samtools' { $snp_workflow = CoGe::Builder::SNP::Samtools::build($snp_params); }
                case 'platypus' { $snp_workflow = CoGe::Builder::SNP::Platypus::build($snp_params); }
                case 'gatk'     { $snp_workflow = CoGe::Builder::SNP::GATK::build($snp_params); }
                else            { die "unknown SNP method"; }
            }
            push @tasks, @{$snp_workflow->{tasks}};
            push @done_files, @{$snp_workflow->{done_files}};
        }
        
        # Add methylation workflow (if specified)
        if ( $self->params->{methylation_params} ) {
            my $method = $self->params->{methylation_params}->{method};
            my $methylation_params = {
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                bam_file => $bam_file,
                raw_bam_file => $raw_bam_files[0], # mdb added 2/29/16 for Bismark, COGE-706
                metadata => $metadata,
                additional_metadata => $additional_metadata,
                read_params => $self->params->{read_params},
                methylation_params => $self->params->{methylation_params},
            };
            
            my $methylation_workflow;
            switch ($method) { #FIXME pass into MeasureMethylation instead
                case 'bismark' { $methylation_workflow = CoGe::Builder::Methylation::Bismark::build($methylation_params); }
                case 'bwameth' { $methylation_workflow = CoGe::Builder::Methylation::BWAmeth::build($methylation_params); }
                else           { die "unknown methylation method"; }
            }
            push @tasks, @{$methylation_workflow->{tasks}};
            push @done_files, @{$methylation_workflow->{done_files}};
            
            # Add metaplot workflow (if specified and genome is annotated)
            if ( $self->params->{methylation_params}->{metaplot_params} ) {
                my $metaplot_workflow = CoGe::Builder::Methylation::Metaplot::build($methylation_params);
                push @tasks, @{$metaplot_workflow->{tasks}};
                push @done_files, @{$metaplot_workflow->{done_files}};
            }
        }
        
        # Add ChIP-seq workflow (if specified)
        if ( $self->params->{chipseq_params} ) {
            my $chipseq_params = {
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_files => $bam_files,
                metadata => $metadata,
                additional_metadata => $additional_metadata,
                read_params => $self->params->{read_params},
                chipseq_params => $self->params->{chipseq_params},
            };
            
            my $chipseq_workflow = CoGe::Builder::Protein::ChIPseq::build($chipseq_params);
            push @tasks, @{$chipseq_workflow->{tasks}};
            push @done_files, @{$chipseq_workflow->{done_files}};
        }
        
        # Add ChIP-seq workflow (if specified)
        if ( $self->params->{chipseq_params} ) {
            my $chipseq_params = {
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_files => $bam_files,
                metadata => $metadata,
                additional_metadata => $additional_metadata,
                read_params => $self->params->{read_params},
                chipseq_params => $self->params->{chipseq_params},
            };
            
            my $chipseq_workflow = CoGe::Builder::Protein::ChIPseq::build($chipseq_params);
            push @tasks, @{$chipseq_workflow->{tasks}};
            push @done_files, @{$chipseq_workflow->{done_files}};
        }
    }
    # Else, all other file types
    else {
        # Generate additional metadata for resulting experiments
        my $annotations = CoGe::Core::Metadata::to_annotations($additional_metadata);
        
        my $input_file = $input_files[0];
        
        # Add conversion step for BigWig files
        if ( $file_type eq 'bw' ) {
            my $wig_task = create_bigwig_to_wig_job(
                staging_dir => $self->staging_dir,
                input_file => $input_file
            );
            push @tasks, $wig_task;
            push @done_files, $wig_task->{outputs}->[1];
            $input_file = $wig_task->{outputs}->[0];  
        }
        
        # Submit workflow to generate experiment
        my $load_task = create_load_experiment_job(
            user => $self->user,
            staging_dir => $self->staging_dir,
            wid => $self->workflow->id,
            gid => $genome->id,
            input_file => $input_file,
            metadata => $metadata,
            annotations => $annotations,
            normalize => $self->params->{normalize} ? $self->params->{normalize_method} : 0
        );
        push @tasks, $load_task;
        push @done_files, $load_task->{outputs}->[1];
    }
    
    # Create notebook
    if ($self->params->{notebook} || $self->params->{notebook_id}) {
        #TODO add_items_to_notebook_job and create_notebook_job and their respective scripts can be consolidated
        if ($self->params->{notebook_id}) { # use existing notebook
            my $t = add_items_to_notebook_job(
                notebook_id => $self->params->{notebook_id},
                user => $self->user, 
                wid => $self->workflow->id,
                staging_dir => $self->staging_dir,
                done_files => \@done_files
            );
            push @tasks, $t;
            push @done_files, $t->{outputs}->[1];
        }
        else { # create new notebook
            my $t = create_notebook_job(
                user => $self->user,
                wid => $self->workflow->id,
                metadata => $metadata,
                staging_dir => $self->staging_dir,
                done_files => \@done_files
            );
            push @tasks, $t;
            push @done_files, $t->{outputs}->[1];
        }
    }
    
    # Send notification email #TODO move into shared module
	if ( $self->params->{email} ) {
	    # Build message body
	    my $body = 'Experiment "' . $metadata->{name} . '" has finished loading.';
        $body .= "\nLink: " . $self->site_url if $self->site_url;
        $body .= "\n\nNote: you received this email because you submitted a job on " .
            "CoGe (http://genomevolution.org) and selected the option to be emailed " .
            "when finished.";
	    
		# Create task
		push @tasks, send_email_job(
			to => $self->user->email,
			subject => 'CoGe Load Experiment done',
			body => $body,
			staging_dir => $self->staging_dir,
			done_files => \@done_files
		);
	}

#    print STDERR Dumper \@tasks, "\n";
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

1;
