package CoGe::Builder::Common::Alignment;

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Accessory::Utils qw(is_fastq_file to_filename detect_paired_end);
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
    my $input_files = $opts{input_files}; # array of file paths
    my $genome      = $opts{genome};
    my $metadata    = $opts{metadata};
    my $load_id     = $opts{load_id};
    my $alignment_params = $opts{alignment_params};
    my $trimming_params  = $opts{trimming_params};
    
    my @tasks;
    
    # Setup paths to data files
    #FIXME this is for LoadExperiment, also need to handle IRODS/FTP/HTTP data from API
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);

    # Check multiple files (if more than one file then all should be FASTQ)
    my $numFastq = 0;
    foreach (@$input_files) {
        $numFastq++ if (is_fastq_file($_));
    }
    if ($numFastq > 0 and $numFastq != @$input_files) {
        my $error = 'Unsupported combination of file types';
        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
        return { error => $error };
    }
    if ($numFastq == 0 and @$input_files > 1) {
        my $error = 'Too many files';
        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
        return { error => $error };
    }
    
    # Decompress and validate the fastq input files
    my (@decompressed, @validated);
    foreach my $input_file (@$input_files) {
        # Decompress
        my $done_file;
        if ( $input_file =~ /\.gz$/ ) {
            push @tasks, create_gunzip_job($input_file);
            $input_file =~ s/\.gz$//;
            $done_file = "$input_file.decompressed";
        }
        push @decompressed, $input_file;
        
        # Validate
        my $validate_task = create_validate_fastq_job($input_file, $done_file);
        push @validated, @{$validate_task->{outputs}}[0];
        push @tasks, $validate_task;
    }
        
    # Trim the fastq input files
    my @trimmed;
    my $trim_reads = 1; #TODO add this as an option in LoadExperiment interface
    if ($trim_reads && $trimming_params) {
        if ($alignment_params->{read_type} eq 'paired') { # mdb added 5/8/15 COGE-624 - enable paired-end support in cutadapt
            # Separate files based on last occurrence of _R1 or _R2 in filename
            my ($m1, $m2) = detect_paired_end($input_files);
            unless (@$m1 and @$m2 and @$m1 == @$m2) {
                my $error = 'Mispaired FASTQ files, m1=' . @$m1 . ' m2=' . @$m2;
                print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
                print STDERR 'm1: ', join(' ', @$m1), "\n", 'm2: ', join(' ', @$m2), "\n";
                return { error => $error };
            }
            
            # Create cutadapt task for each file pair
            for (my $i = 0;  $i < @$m1;  $i++) { 
                my $file1 = shift @$m1;
                my $file2 = shift @$m2;
                my $trim_task = create_cutadapt_job(
                    fastq => [ $file1, $file2 ],
                    validated => [ "$file1.validated", "$file2.validated" ],
                    staging_dir => $staging_dir,
                    params => $trimming_params
                );
                push @trimmed, @{$trim_task->{outputs}};
                push @tasks, $trim_task;
            }
        }
        else { # single-ended
            # Create cutadapt task for each file
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
    else { # no trimming
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
        my $error = "Unrecognized alignment tool '" . $alignment_params->{tool};
        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
        return { error => $error };
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

    # Get custom metadata to add to experiment
    my $additional_md = generate_additional_metadata($trimming_params, $alignment_params);

    # Load alignment
    my $load_task = create_load_bam_job(
        user => $user,
        metadata => $metadata,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        annotations => $additional_md,
        wid => $wid,
        gid => $gid,
        bam_file => $sorted_bam_file
    );
    push @tasks, $load_task;
    
    return {
        tasks => \@tasks,
        bam_file => $sorted_bam_file,
        done_files => [
            $sorted_bam_file,
            $load_task->{outputs}->[1]
        ]
    }
}

sub generate_additional_metadata {
    my ($trimming_params, $alignment_params) = @_;
    my @annotations;
    
    push @annotations, qq{https://genomevolution.org/wiki/index.php/Expression_Analysis_Pipeline||note|Generated by CoGe's RNAseq Analysis Pipeline};
    
    if ($trimming_params) {
        push @annotations, 'note|cutadapt '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--quality-base', '-m'));
    }

    if ($alignment_params && $alignment_params->{tool}) {
        if ($alignment_params->{tool} eq 'tophat') { # tophat
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|tophat ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-g'));
        }
        else { # gsnap
            push @annotations, qq{note|gmap_build};
            push @annotations, 'note|gsnap ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-n', '-Q', '--gap-mode', '--nofails'));
        }
    }
    
    return \@annotations;
}

1;