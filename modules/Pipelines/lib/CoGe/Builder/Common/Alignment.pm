package CoGe::Builder::Common::Alignment;

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Accessory::Utils qw(is_fastq_file to_filename);
use CoGe::Builder::CommonTasks;

our $CONF = CoGe::Accessory::Web::get_defaults();

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT_OK = qw(build);
}

sub build {
    my %opts = @_;
    my $user        = $opts{user};
    my $wid         = $opts{wid};
    my $input_files = $opts{input_files};
    my $genome      = $opts{genome};
    my $metadata    = $opts{metadata};
    my $additional_metadata = $opts{additional_metadata};
    my $options     = $opts{options};
    my $alignment_params = $opts{alignment_params};
    my $trimming_params  = $opts{trimming_params};
    
    my @tasks;
    
    # Setup paths to data files
    #FIXME this is for LoadExperiment, also need to handle IRODS/FTP/HTTP data from API
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $upload_dir = get_upload_path($user->name, $options->{load_id});
    my @files = map { catfile($upload_dir, $_->{path}) } @$input_files;

    # Check multiple files (if more than one file then all should be FASTQ)
    my $numFastq = 0;
    foreach (@files) {
        $numFastq++ if (is_fastq_file($_));
    }
    return encode_json({ error => 'Unsupported combination of file types' }) if ($numFastq > 0 and $numFastq != @files);
    return encode_json({ error => 'Too many files' }) if ($numFastq == 0 and @files > 1);
    
    # Decompress the fastq input files (when necessary)
    my @decompressed;
    foreach my $file (@files) {
        my $done_file;
        if ( $file =~ /\.gz$/ ) {
            push @tasks, create_gunzip_job($file);
            $file =~ s/\.gz$//;
        }
        push @decompressed, $file;
    }
        
    # Validate the fastq input files
    my @validated;
    foreach my $file (@decompressed) {
        my $validate_task = create_validate_fastq_job($file, "$file.decompressed");
        push @validated, @{$validate_task->{outputs}}[0];
        push @tasks, $validate_task;
    }
    
    # Trim the fastq input files
    my @trimmed;
    my $trim_reads = 1; #TODO add this as an option in LoadExperiment interface
    if ($trim_reads) {
        if ($alignment_params->{read_type} eq 'paired') { # mdb added 5/8/15 COGE-624 - enable paired-end support in cutadapt
            my @m1 = sort grep { $_ =~ /\_R1/ } @files;
            my @m2 = sort grep { $_ =~ /\_R2/ } @files;
            unless (@m1 and @m2 and @m1 == @m2) {
                return encode_json({ error => 'Mis-paired FASTQ files' });
            }
            
            for (my $i = 0;  $i < @m1;  $i++) {
                my $file1 = shift @m1;
                my $file2 = shift @m2;
                my $trim_task = create_cutadapt_job(
                    fastq => [ $file1, $file2 ],
                    validated => [ "$file1.validated", "$file2.validated" ]
                    staging_dir => $staging_dir,
                    params => $trimming_params
                );
                push @trimmed, $trim_task->{outputs}->[0];
                push @tasks, $trim_task;
            }
        }
        else { # assume single-ended
            foreach my $file (@decompressed) {
                my $trim_task = create_cutadapt_job(
                    fastq => $file,
                    validated => "$file.validated",
                    staging_dir => $staging_dir,
                    params => $trimming_params
                );
                push @trimmed, $trim_task->{outputs}->[0];
                push @tasks, $trim_task;
            }
        }
    }
    else {
        push @trimmed, @decompressed;
    }

    # Get genome cache path
    my $gid = $genome->id;
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
    
    # Add aligner workflow
    my ($alignment_tasks, $alignment_results);
    if ($alignment_params->{tool} eq 'tophat') {
        # Generate gff if genome annotated
        my $gff_file;
        if ( $genome->has_gene_features ) {
            my $gff_task = create_gff_generation_job(
                gid => $gid,
                organism_name => $genome->organism->name,
                #validated => "$fastq.validated"
            );
            $gff_file = $gff_task->{outputs}->[0];
            push @tasks, $gff_task;
        }
        
        ($alignment_tasks, $alignment_results) = create_tophat_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $alignment_params->{read_type},
            gff => $gff_file,
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    elsif ($alignment_params->{tool} eq 'gsnap') {
        ($alignment_tasks, $alignment_results) = create_gsnap_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $alignment_params->{read_type},
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    else {
        print STDERR "Error: unrecognized alignment tool '", $alignment_params->{tool}, "'\n";
        return;
    }
    push @tasks, @$alignment_tasks;

    # Sort and index the bam output file
    my $bam_file = $alignment_results->{bam_file};
    my $sort_bam_task = create_bam_sort_job(
        input_file => $bam_file, 
        staging_dir => $staging_dir
    );
    push @tasks, $sort_bam_task;
    my $sorted_bam_file = $sort_bam_task->{outputs}->[0];
    
    push @tasks, create_bam_index_job(
        input_file => $sorted_bam_file
    );

    # Load alignment
    my $load_task = create_load_bam_job(
        user => $user,
        metadata => $metadata,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        annotations => $additional_metadata,
        wid => $wid,
        gid => $gid,
        bam_file => $sorted_bam_file
    );
    push @tasks, $load_task;
    
    # Save outputs for retrieval by downstream tasks
    my %results = (
        bam_file => $sorted_bam_file,
        metadata => generate_additional_metadata($trimming_params, $alignment_params),
        done_files => [
            $sorted_bam_file,
            $load_task->{outputs}->[1]
        ]
    );
    
    return (\@tasks, \%results);
}

sub generate_additional_metadata {
    my ($trimming_params, $alignment_params) = @_;
    my @annotations;
    
    if ($trimming_params) {
        push @annotations, 'note|cutadapt '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--quality-base', '-m'));
    }

    if ($alignment_params->{tool}) {
        if ($alignment_params->{tool} eq 'tophat') { # tophat
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|tophat ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-g'));
        }
        else { # gsnap
            push @annotations, qq{note|gmap_build};
            push @annotations, 'note|gsnap ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-n', '-Q', '--gap-mode', '--nofails'));
        }
    }
    
    return \@annotations;
}

1;