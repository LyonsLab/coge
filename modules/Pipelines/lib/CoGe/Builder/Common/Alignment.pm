package CoGe::Builder::Common::Alignment;

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Core::Metadata qw(to_annotations);
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
    my $additional_metadata = $opts{additional_metadata};
    my $load_id     = $opts{load_id};
    my $params      = $opts{params};
    my $read_params      = $params->{read_params};
    my $trimming_params  = $params->{trimming_params};
    my $alignment_params = $params->{alignment_params};
    my $chipseq_params   = $params->{chipseq_params};
    
    my @tasks;
    
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);

# mdb removed 11/6/15 COGE-673
#    # Check multiple files (if more than one file then all should be FASTQ)
#    my $numFastq = 0;
#    foreach (@$input_files) {
#        $numFastq++ if (is_fastq_file($_));
#    }
#    if ($numFastq > 0 and $numFastq != @$input_files) {
#        my $error = 'Unsupported combination of file types';
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }
#    if ($numFastq == 0 and @$input_files > 1) {
#        my $error = 'Too many files';
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }
    
    # Decompress and validate the fastq input files
    my (@decompressed, @validated);
    foreach my $input_file (@$input_files) {
        # Decompress
        my $done_file = "$input_file.done";
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
    if ($trimming_params) {
        my ($fastq1, $fastq2);
        if ($read_params->{read_type} eq 'paired') {
            # Separate files based on last occurrence of _R1 or _R2 in filename
            my ($m1, $m2) = detect_paired_end(\@decompressed);
            unless (@$m1 and @$m2 and @$m1 == @$m2) {
                my $error = 'Mispaired FASTQ files, m1=' . @$m1 . ' m2=' . @$m2;
                print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
                print STDERR 'm1: ', join(' ', @$m1), "\n", 'm2: ', join(' ', @$m2), "\n";
                return { error => $error };
            }
            $fastq1 = $m1;
            $fastq2 = $m2;
        }
        else { # default to single-ended
            $fastq1 = \@decompressed;
        }
        
        my %params = (
            fastq1 => $fastq1,
            fastq2 => $fastq2,
            validated => \@validated,
            staging_dir => $staging_dir,
            read_params => $read_params,
            trimming_params => $trimming_params        
        );
        
        my ($tasks, $outputs);
        if ($trimming_params->{trimmer} eq 'cutadapt') {
            ($tasks, $outputs) = create_cutadapt_workflow(%params);
        }
        elsif ($trimming_params->{trimmer} eq 'trimgalore') {
            ($tasks, $outputs) = create_trimgalore_workflow(%params);
        }
        push @trimmed, @$outputs;
        push @tasks, @$tasks;
    }
    else { # no trimming
        push @trimmed, @decompressed;
    }

    # Get genome cache path
    my $gid = $genome->id;
    my $fasta_cache_dir = get_genome_cache_path($gid);

    # Reheader the fasta file
    my $fasta = get_genome_file($gid);
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta,
        cache_dir => $fasta_cache_dir
    );
    my $reheader_fasta = $tasks[-1]->{outputs}[0];

    # Index the fasta file
    push @tasks, create_fasta_index_job(
        fasta => $reheader_fasta
    );
    
    my %params = ( 
        fasta => catfile($fasta_cache_dir, $reheader_fasta),
        fastq => \@trimmed,
        validated => \@validated,
        gid => $gid,
        encoding => $read_params->{encoding},
        read_type => $read_params->{read_type},
        staging_dir => $staging_dir,
        params => $alignment_params,
    );
    
    # Add aligner workflow
    my ($alignment_tasks, $alignment_results);
    $alignment_params = {} unless $alignment_params;
    
    if ($alignment_params->{tool} eq 'hisat2') {
        ($alignment_tasks, $alignment_results) = create_hisat2_workflow(%params);
    }
    elsif ($alignment_params->{tool} eq 'bowtie2') {
        if ($chipseq_params) {
            ($alignment_tasks, $alignment_results) = create_bowtie2_chipseq_workflow(%params);
        }
        else {
            ($alignment_tasks, $alignment_results) = create_bowtie2_workflow(%params);
        }
    }
    elsif ($alignment_params->{tool} eq 'tophat') {
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
        
        $params{gff} = $gff_file;
        ($alignment_tasks, $alignment_results) = create_tophat_workflow(%params);
    }
    elsif ($alignment_params->{tool} eq 'bismark') {
        ($alignment_tasks, $alignment_results) = create_bismark_workflow(%params);
    }
    elsif ($alignment_params->{tool} eq 'bwameth') {
        ($alignment_tasks, $alignment_results) = create_bwameth_workflow(%params);
    }
    else { # ($alignment_params->{tool} eq 'gsnap') { # default
        ($alignment_tasks, $alignment_results) = create_gsnap_workflow(%params);
    }
# mdb removed 9/15/15 -- make GSNAP the default aligner
#    else {
#        my $error = "Unrecognized alignment tool '" . $alignment_params->{tool};
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }
    
    unless (@$alignment_tasks && $alignment_results) {
        print STDERR "CoGe::Builder::Common::Alignment ERROR: Invalid alignment workflow\n";
        return { error => 'Invalid alignment workflow' };
    }
    push @tasks, @$alignment_tasks;

    my @bam_files;
    push @bam_files, $alignment_results->{bam_file} if ($alignment_results->{bam_file});
    push @bam_files, @{$alignment_results->{bam_files}} if ($alignment_results->{bam_files});
    
    my (@sorted_bam_files, @load_task_outputs);
    foreach my $bam_file (@bam_files) {
        # Sort and index the bam output file(s)
        my $sort_bam_task = create_bam_sort_job(
            input_file => $bam_file, 
            staging_dir => $staging_dir
        );
        push @tasks, $sort_bam_task;
        my $sorted_bam_file = $sort_bam_task->{outputs}->[0];
        push @sorted_bam_files, $sorted_bam_file;
        
        push @tasks, create_bam_index_job(
            input_file => $sorted_bam_file
        );
    
        # Get custom metadata to add to experiment
        my $annotations = generate_additional_metadata($read_params, $trimming_params, $alignment_params);
        my @annotations2 = CoGe::Core::Metadata::to_annotations($additional_metadata);
        push @$annotations, @annotations2;
    
        # Load alignment
        my $load_task = create_load_bam_job(
            user => $user,
            metadata => $metadata,
            staging_dir => $staging_dir,
            result_dir => $result_dir,
            annotations => $annotations,
            wid => $wid,
            gid => $gid,
            bam_file => $sorted_bam_file
        );
        push @tasks, $load_task;
        push @load_task_outputs, $load_task->{outputs}->[1];
    }
    
    return {
        tasks => \@tasks,
        bam_files => \@sorted_bam_files,
        done_files => [
            @sorted_bam_files,
            @load_task_outputs
        ]
    }
}

sub generate_additional_metadata {
    my ($read_params, $trimming_params, $alignment_params) = @_;
    my @annotations;
    
    push @annotations, qq{https://genomevolution.org/wiki/index.php/Expression_Analysis_Pipeline||note|Generated by CoGe's RNAseq Analysis Pipeline};
    
    if ($trimming_params && $trimming_params->{trimmer}) {
        if ($trimming_params->{trimmer} eq 'cutadapt') {
            push @annotations, 'note|cutadapt '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--quality-base', '-m'));
        }
        elsif ($trimming_params->{trimmer} eq 'trimgalore') {
            push @annotations, 'note|trimgalore '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--length', '-a'));
        }
    }

    if ($alignment_params && $alignment_params->{tool}) {
        if ($alignment_params->{tool} eq 'hisat2') {
            push @annotations, qq{note|hisat2_build};
            my $params = ($read_params->{encoding} eq '64' ? '--phred64' : '--phred33');
            push @annotations, 'note|hisat2 ' . $params;
        }
        elsif ($alignment_params->{tool} eq 'tophat') {
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|tophat ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-g'));
        }
        elsif ($alignment_params->{tool} eq 'bismark') {
            push @annotations, qq{note|bismark_genome_preparation};
            push @annotations, 'note|bismark ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-L'));
        }
        elsif ($alignment_params->{tool} eq 'bwameth') {
            push @annotations, qq{note|bwameth index};
            push @annotations, 'note|bwameth';
        }
        else { # gsnap
            push @annotations, qq{note|gmap_build};
            push @annotations, 'note|gsnap ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-n', '-Q', '--gap-mode', '--nofails'));
        }
    }
    
    return \@annotations;
}

1;
