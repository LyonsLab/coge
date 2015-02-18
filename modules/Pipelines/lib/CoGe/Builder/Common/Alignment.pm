package CoGe::Builder::Common::Alignment;

use Moose;

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Accessory::Utils qw(is_fastq_file to_filename);
use CoGe::Builder::CommonTasks;

sub build {
    my $self = shift;
    
    # Validate inputs and set defaults
    unless ($self->options && $self->params) {
        print STDERR "Alignment::build missing options or params\n";
        return;
    }
    
    my $gid = $self->params->{gid};
    unless ($gid) {
        print STDERR "Alignment::build missing gid\n";
        return;
    }
    
    my $data = $self->options->{source_data};
    unless (defined $data && @$data) {
        print STDERR "Alignment::build missing source_data\n";
        return;
    }
    
    my $alignment_params = $self->params->{alignment_params};
    unless ($alignment_params) {
        print STDERR "Alignment::build missing alignment_params\n";
        return;
    }
    
    my $read_type = $self->params->{alignment_params}{read_type};
    unless ($read_type) {
        print STDERR "Alignment::build missing read_type\n";
        return;
    }
    
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        print STDERR "Alignment::build missing metadata\n";
        return;
    }
    
    # Initialize workflow
    unless ($self->workflow) {
        $self->workflow( $self->jex->create_workflow(name => "Align reads to reference", init => 1) );
        return unless $self->workflow->id;
    }
    
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    
    # Build workflow steps -----------------------------------------------------
    my @tasks;
    
    # Setup paths to data files
    #FIXME this is for LoadExperiment, also need to handle IRODS/FTP/HTTP data from API
    my $upload_dir = get_upload_path($self->user->name, $self->options->{load_id});
    my @files = map { catfile($upload_dir, $_->{path}) } @$data;

    # Check multiple files (if more than one file then all should be FASTQ)
    my $numFastq = 0;
    foreach (@files) {
        $numFastq++ if (is_fastq_file($_));
    }
    return return encode_json({ error => 'Unsupported combination of file types' }) if ($numFastq > 0 and $numFastq != @files);
    return return encode_json({ error => 'Too many files' }) if ($numFastq == 0 and @files > 1);
    
    # Process the fastq input files
    my (@validated, @trimmed);
    foreach my $file (@files) {
        # Decompress (if necessary)
        my $done_file;
        if ( $file =~ /\.gz$/ ) {
            push @tasks, create_gunzip_job($file);
            $file =~ s/\.gz$//;
            $done_file = "$file.decompressed";
        }
        
        # Validate
        my $validate_task = create_validate_fastq_job($file, $done_file);
        push @validated, @{$validate_task->{outputs}}[0];
        push @tasks, $validate_task;
    
        # Trim
        my $trim_reads = 1; #FIXME hook this as option in LoadExperiment interface
        if ($trim_reads) {
            my $trim_task = create_cutadapt_job(
                fastq => $file, 
                validated => "$file.validated", 
                staging_dir => $staging_dir,
                params => $self->params->{cutadapt_params}
            );
            push @trimmed, @{$trim_task->{outputs}}[0];
            push @tasks, $trim_task;
        }
        else {
            push @trimmed, $file;
        }
    }

    # Get genome cache path
    my $fasta_cache_dir = get_genome_cache_path($gid);

    # Reheader the fasta file
    my $fasta = get_genome_file($gid);
    my $reheader_fasta = to_filename($fasta) . ".reheader.faa";
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta,
        reheader_fasta => $reheader_fasta,
        cache_dir => $fasta_cache_dir
    );

    # Index the fasta file
    push @tasks, create_fasta_index_job(
        fasta => catfile($fasta_cache_dir, $reheader_fasta),
        cache_dir => $fasta_cache_dir
    );
    
    # Get genome (needed for organism name next)
    my $genome = $self->db->resultset('Genome')->find($gid);
    return unless ($genome);
    # TODO add permissions check here -- or will it happen in Request::Genome?

    # Generate gff if genome annotated
    my $gff_file;
    if ( $genome->has_gene_features ) {
        my $gff_task = create_gff_generation_job(
            gid => $gid,
            organism_name => $genome->organism->name,
            #validated => "$fastq.validated"
        );
        $gff_file = @{$gff_task->{outputs}}[0];
        push @tasks, $gff_task;
    }

    # Add aligner workflow
    my ($bam, @alignment_tasks);
    if ($alignment_params->{tool} eq 'tophat') {
        ($bam, @alignment_tasks) = create_tophat_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $read_type,
            gff => $gff_file,
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    elsif ($alignment_params->{tool} eq 'gsnap') {
        ($bam, @alignment_tasks) = create_gsnap_workflow(
            gid => $gid,
            fasta => $reheader_fasta,
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $read_type,
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    else {
        print STDERR "Error: unrecognized alignment tool '", $alignment_params->{tool}, "'\n";
        return;
    }
    push @tasks, @alignment_tasks;

    # Load alignment
    push @tasks, create_load_bam_job(
        user => $self->user,
        metadata => $metadata,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        annotations => '', #FIXME 12/12/14
        wid => $self->workflow->id,
        gid => $gid,
        bam_file => $bam
    );
    
    # Add all tasks to workflow
    $self->workflow->add_jobs(\@tasks);

    # Save outputs for retrieval by downstream tasks
    $self->outputs->{reheader_fasta} = $reheader_fasta;
    $self->outputs->{bam} = $bam;
    $self->outputs->{gff} = $gff_file;
    #return (\@tasks, $bam, $reheader_fasta, $gff_file);
    return 1;
}

with qw(CoGe::Builder::Buildable);

1;