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
    my $read_params      = $opts{read_params};
    my $trimming_params  = $opts{trimming_params};
    my $alignment_params = $opts{alignment_params};
    
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
        
        if ($trimming_params->{trimmer} eq 'cutadapt') {
            my ($tasks, $outputs) = create_cutadapt_workflow(
                fastq1 => $fastq1,
                fastq2 => $fastq2,
                validated => \@validated,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @trimmed, @$outputs;
            push @tasks, @$tasks;
        }
        elsif ($trimming_params->{trimmer} eq 'trimgalore') {
            my ($tasks, $outputs) = create_trimgalore_workflow(
                fastq1 => $fastq1,
                fastq2 => $fastq2,
                validated => \@validated,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @trimmed, @$outputs;
            push @tasks, @$tasks;
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
    if ($alignment_params && $alignment_params->{tool} eq 'hisat2') {
        ($alignment_tasks, $alignment_results) = create_hisat2_workflow(
            fasta => $fasta,
        	fastq => \@trimmed,
            gid => $gid,
            encoding => $read_params->{encoding},
            read_type => $read_params->{read_type},
            staging_dir => $staging_dir
        );    	
    }
    elsif ($alignment_params && $alignment_params->{tool} eq 'tophat') {
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
            encoding => $read_params->{encoding},
            read_type => $read_params->{read_type},
            gff => $gff_file,
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    elsif ($alignment_params && $alignment_params->{tool} eq 'bismark') {
        ($alignment_tasks, $alignment_results) = create_bismark_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            encoding => $read_params->{encoding},
            read_type => $read_params->{read_type},
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    elsif ($alignment_params && $alignment_params->{tool} eq 'bwameth') {
        ($alignment_tasks, $alignment_results) = create_bwameth_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            #encoding => $read_params->{encoding}, # bwameth doesn't have encoding option, must auto-detect?
            read_type => $read_params->{read_type},
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    else { # ($alignment_params->{tool} eq 'gsnap') { # default
        ($alignment_tasks, $alignment_results) = create_gsnap_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $read_params->{read_type},
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
# mdb removed 9/15/15 -- make GSNAP the default aligner
#    else {
#        my $error = "Unrecognized alignment tool '" . $alignment_params->{tool};
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }
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
            my $params = join(' ', map { $_.' '.$alignment_params->{$_} } ('-p', '-x', '-S'));
            if ($read_params->{encoding} eq '64') {
            	$params .= ' --phred64';
            }
            else {
            	$params .= ' --phred33';
            }
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
