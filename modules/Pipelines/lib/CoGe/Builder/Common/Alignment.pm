package CoGe::Builder::Common::Alignment;

use strict;
use warnings;

use Switch;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use Clone qw(clone);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Core::Metadata qw(to_annotations);
use CoGe::Accessory::Utils qw(is_fastq_file to_filename detect_paired_end to_filename_base to_filename_without_extension);
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
    my $input_files = $opts{input_files}; # array of paths of FASTQ files
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
    my (@trimmed, @trimming_done_files);
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
            done_files => \@validated,
            staging_dir => $staging_dir,
            read_params => $read_params,
            trimming_params => $trimming_params
        );
        
        my ($tasks, $outputs, $done_files);
        if ($trimming_params->{trimmer} eq 'cutadapt') {
            ($tasks, $outputs, $done_files) = create_cutadapt_workflow(%params);
        }
        elsif ($trimming_params->{trimmer} eq 'trimgalore') {
            ($tasks, $outputs, $done_files) = create_trimgalore_workflow(%params);
        }
        push @trimmed, @$outputs;
        push @tasks, @$tasks;
        push @trimming_done_files, @$done_files if $done_files; # kludge for JEX
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
    
    # Add aligner workflow
    my %params = ( 
        done_files => [ @validated, @trimming_done_files ],
        fasta => $reheader_fasta,
        fastq => \@trimmed,
        gid => $gid,
        encoding => $read_params->{encoding},
        read_type => $read_params->{read_type},
        staging_dir => $staging_dir,
        params => $alignment_params
    );
    $params{doSeparately} = 1 if $chipseq_params;
    
    my ($alignment_tasks, $alignment_results);
    $alignment_params = {} unless $alignment_params;

    switch( _get_aligner($alignment_params) ) {
        case 'hisat2'  { ($alignment_tasks, $alignment_results) = create_hisat2_workflow(%params) }
        case 'bowtie2' { ($alignment_tasks, $alignment_results) = create_bowtie2_workflow(%params) }
        case 'tophat'  {
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
        case 'bismark' { ($alignment_tasks, $alignment_results) = create_bismark_workflow(%params) }
        case 'bwameth' { ($alignment_tasks, $alignment_results) = create_bwameth_workflow(%params) }
        case 'bwa'     { ($alignment_tasks, $alignment_results) = create_bwa_workflow(%params) }
        case 'gsnap'   { ($alignment_tasks, $alignment_results) = create_gsnap_workflow(%params) }
    }

    unless (@$alignment_tasks && $alignment_results) {
        print STDERR "CoGe::Builder::Common::Alignment ERROR: Invalid alignment workflow\n";
        return { error => 'Invalid alignment workflow' };
    }
    push @tasks, @$alignment_tasks;

    my @raw_bam_files;
    push @raw_bam_files, $alignment_results->{bam_file} if ($alignment_results->{bam_file});
    push @raw_bam_files, @{$alignment_results->{bam_files}} if ($alignment_results->{bam_files});
    
    my (@sorted_bam_files, @load_task_outputs);
    foreach my $bam_file (@raw_bam_files) {
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
        
        # Add bam filename to experiment name for ChIP-seq pipeline
        my $md = clone($metadata);
        if (@raw_bam_files > 1) {
            $md->{name} .= ' (' . to_filename_base($sorted_bam_file) . ')';
        }
    
        # Load alignment
        my $load_task = create_load_bam_job(
            user => $user,
            metadata => $md,
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
        raw_bam_files => \@raw_bam_files, # mdb added 2/29/16 for Bismark, COGE-706
        done_files => [
            @sorted_bam_files,
            @load_task_outputs
        ]
    }
}

sub create_validate_fastq_job {
    my $fastq = shift;
    my $done_file = shift;

    my $cmd = catfile($CONF->{SCRIPTDIR}, "validate_fastq.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $inputs = [ $fastq ];
    push @$inputs, $done_file if ($done_file);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fastq, 0] # mdb changed 3/1/16 from 1 to 0, COGE-707
        ],
        inputs => $inputs,
        outputs => [
            "$fastq.validated"
        ],
        description => "Validating " . basename($fastq)
    };
}

sub create_cutadapt_job {
    my %opts = @_;

    # Required params
    my $fastq = $opts{fastq};            # for single fastq file (backwards compatibility) or two paired-end fastq files (new functionality)
    my $done_files = $opts{done_files};  # input dependency from previous task, one or two files based on fastq arg
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $trimming_params = $opts{trimming_params} // {};
    my $q = $trimming_params->{'-q'} // 25;
    my $m = $trimming_params->{'-m'} // 17;
    my $read_params = $opts{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');
    $done_files = [ $done_files ] unless (ref($done_files) eq 'ARRAY');

    my $name = join(', ', map { basename($_) } @$fastq);
    my @inputs = ( @$fastq, @$done_files);
    my @outputs = map { catfile($staging_dir, to_filename($_) . '.trimmed.fastq') } @$fastq;
    my @done_files = map { $_ . '.done' } @outputs;

    # Build up command/arguments string
    my $cmd = get_command_path('CUTADAPT');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $arg_str;
    $arg_str .= $cmd . ' ';
    $arg_str .= "-q $q --quality-base=$encoding -m $m -o $outputs[0] ";
    $arg_str .= "-p $outputs[1] " if ($read_type eq 'paired'); # paired-end

    return {
        cmd => catfile($CONF->{SCRIPTDIR}, 'cutadapt.pl'), # this script was created because JEX can't handle Cutadapt's paired-end argument syntax
        script => undef,
        args => [
            [$read_type, '', 0],
            [$staging_dir, '', 0],
            ['"'.$arg_str.'"', '', 0],
            ['', join(' ', @$fastq), 0]
        ],
        inputs => \@inputs,
        outputs => \@outputs,
        done_files => \@done_files, # JEX will ignore this
        description => "Trimming (cutadapt) $name"
    };
}

sub create_cutadapt_workflow {
    my %opts = @_;
    my $fastq1 = $opts{fastq1}; # array ref of left reads (or all reads if single-ended)
    my $fastq2 = $opts{fastq2}; # array ref of right reads (or undef if single-ended)
    my $done_files = $opts{done_files};
    my $staging_dir = $opts{staging_dir};
    my $read_params = $opts{read_params};
    my $read_type = $read_params->{read_type} // 'single'; #/
    my $trimming_params = $opts{trimming_params};

    my (@tasks, @outputs, @done_files);

    if ($read_type eq 'single') { # single-ended
        # Create cutadapt task for each file
        foreach my $file (@$fastq1) {
            my $task = create_cutadapt_job(
                fastq => $file,
                done_files => $done_files,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, $task->{outputs}->[0];
            push @done_files, @{$task->{done_files}};
            push @tasks, $task;
        }
    }
    else { # paired-end
        # Create cutadapt task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my $file1 = shift @$fastq1;
            my $file2 = shift @$fastq2;
            my $task = create_cutadapt_job(
                fastq => [ $file1, $file2 ],
                done_files => $done_files,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, @{$task->{outputs}};
            push @done_files, @{$task->{done_files}};
            push @tasks, $task;
        }
    }

    return ( \@tasks, \@outputs, \@done_files );
}

sub create_trimgalore_job {
    my %opts = @_;

    # Required params
    my $fastq = $opts{fastq};            # for single fastq file (backwards compatibility) or two paired-end fastq files (new functionality)
    my $done_files = $opts{done_files};  # input dependency from previous task, one or two files based on fastq arg
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $trimming_params = $opts{trimming_params} // {}; #/
    my $q = $trimming_params->{'-q'} // 20; #/
    my $length = $trimming_params->{'-length'} // 20; #/
    my $a = $trimming_params->{'-a'};
    my $read_params = $opts{read_params} // {}; #/
    my $encoding = $read_params->{encoding} // 33; #/
    my $read_type = $read_params->{read_type} // 'single'; #/

    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');
    $done_files = [ $done_files ] unless (ref($done_files) eq 'ARRAY');

    my $name = join(', ', map { basename($_) } @$fastq);
    my @inputs = ( @$fastq, @$done_files);

    # Build up command/arguments string
    my $cmd = $CONF->{TRIMGALORE} || 'trim_galore';
    $cmd = 'nice ' . $cmd; # run at lower priority

    # Create staging dir
    $cmd = "mkdir -p $staging_dir && " . $cmd;

    my $args = [
        ['--output_dir', $staging_dir, 0],
        ['-q', $q, 0],
        ['--length', $length, 0],
    ];

    my $phred = ($encoding == 64 ? '--phred64' : '--phred33');
    push @$args, [$phred, '', 0];

    if (defined $a) {
        push @$args, ['-a', '"'.$a.'"', 0];
    }

    my @outputs;
    if ($read_type eq 'paired') {
        push @$args, ['--paired', join(' ', @$fastq), 0];

        my ($r1, $r2) = @$fastq;
        @outputs = ( catfile($staging_dir, to_filename_without_extension($r1) . '_val_1.fq'),
                     catfile($staging_dir, to_filename_without_extension($r2) . '_val_2.fq') );
    }
    else { # single
        my ($file) = @$fastq;
        push @$args, ['', $file, 0];
        @outputs = ( catfile($staging_dir, to_filename_without_extension($file) . '_trimmed.fq') );
    }

    # kludge to fix JEX sequencing
    foreach (@$args) {
        $cmd .= ' ' . $_->[0] . ' ' . $_->[1];
    }
    my @done_files = map { $_ . '.done' } @outputs;
    foreach (@done_files) {
        $cmd .= " && touch $_";
    }

    return {
        cmd => $cmd,
        script => undef,
        args => [], #$args
        inputs => \@inputs,
        outputs => \@outputs,
        done_files => \@done_files,
        description => "Trimming (trimgalore) $name"
    };
}

sub create_trimgalore_workflow {
    my %opts = @_;
    my $fastq1 = $opts{fastq1} // []; # array ref of left reads (or all reads if single-ended)
    my $fastq2 = $opts{fastq2} // []; # array ref of right reads (or undef if single-ended)
    my $done_files = $opts{done_files};
    my $staging_dir = $opts{staging_dir};
    my $read_params = $opts{read_params};
    my $trimming_params = $opts{trimming_params};

    my (@tasks, @outputs, @done_files);

    if (!$read_params->{read_type} || $read_params->{read_type} eq 'single') { # single-ended
        # Create trimgalore task for each file
        foreach my $file (@$fastq1, @$fastq2) {
            my $task = create_trimgalore_job(
                fastq => $file,
                done_files => $done_files,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, $task->{outputs}->[0];
            push @done_files, @{$task->{done_files}};
            push @tasks, $task;
        }
    }
    else { # paired-end
        # Create trimgalore task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my $file1 = shift @$fastq1;
            my $file2 = shift @$fastq2;
            my $task = create_trimgalore_job(
                fastq => [ $file1, $file2 ],
                done_files => $done_files,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, @{$task->{outputs}};
            push @done_files, @{$task->{done_files}};
            push @tasks, $task;
        }
    }

    return ( \@tasks, \@outputs, \@done_files );
}

sub create_hisat2_workflow {
    my %opts = @_;
    my $fasta        = $opts{fasta};
    my $fastq        = $opts{fastq};
    my $done_files   = $opts{done_files};
    my $gid          = $opts{gid};
    my $encoding     = $opts{encoding};
    my $read_type    = $opts{read_type} // 'single';
    my $staging_dir  = $opts{staging_dir};
    my $doSeparately = $opts{doSeparately}; # for ChIP-seq pipeline, align each fastq in separate runs rather than together

    my @tasks;

    # Add index task
    my $index_task = create_hisat2_index_job($gid, $fasta);
    push @tasks, $index_task;

    # Add one or more alignment tasks
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            my $task = create_hisat2_alignment_job(
                fastq => [ $file ],
                done_files => $done_files,
                gid => $gid,
                index_files => $index_task->{outputs},
                encoding => $encoding,
                read_type => $read_type,
                staging_dir => $staging_dir
            );
            push @tasks, $task;
            push @sam_files, $task->{outputs}->[0];
        }
    }
    else { # standard Bowtie run (all fastq's at once)
        my $task = create_hisat2_alignment_job(
            fastq => $fastq,
            done_files => $done_files,
            gid => $gid,
            index_files => $index_task->{outputs},
            encoding => $encoding,
            read_type => $read_type,
            staging_dir => $staging_dir
        );
        push @tasks, $task;
        push @sam_files, $task->{outputs}->[0];
    }

    # Add one or more sam-to-bam tasks
    my %results;
    foreach my $file (@sam_files) {
        my $task = create_sam_to_bam_job($file, $staging_dir);
        push @tasks, $task;
        push @{$results{bam_files}}, $task->{outputs}->[0];
    }

    return (\@tasks, \%results);
}

sub create_hisat2_alignment_job {
    my %opts = @_;
    my $fastq       = $opts{fastq};
    my $done_files  = $opts{done_files};
    my $gid         = $opts{gid};
    my $index_files = $opts{index_files};
    my $encoding    = $opts{encoding};
    my $read_type   = $opts{read_type};
    my $staging_dir = $opts{staging_dir};

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.sam';

	my $args = [
		['-p', '32', 0],
                ['--dta-cufflinks', '', 0], # mdb added 8/9/16 for Cufflinks error "BAM record error: found spliced alignment without XS attribute"
		['-x', catfile(get_genome_cache_path($gid), 'hisat2_index', 'genome.reheader'), 0],
		['-S', $output_file, 0]
    ];

    if ($encoding eq '64') {
        push $args, ['--phred64', '', 0];
    }
    else {
        push $args, ['--phred33', '', 0];
    }

    if ($read_type eq 'single') {
    	push $args, ['-U', join(',', @$fastq), 0];
    }
    else {
		my ($m1, $m2) = detect_paired_end($fastq);
    	push $args, ['-1', join(',', sort @$m1), 0];
    	push $args, ['-2', join(',', sort @$m2), 0];
	}

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

	return {
        cmd => 'nice ' . get_command_path('HISAT2'),
        script => undef,
        args => $args,
        inputs => [ @$fastq, @$done_files, @$index_files ],
        outputs => [ catfile($staging_dir, $output_file) ],
        description => "Aligning $desc using HISAT2"
	};
}

sub create_hisat2_index_job {
    my $gid = shift;
    my $fasta = shift;

    my $cache_dir = catdir(get_genome_cache_path($gid), "hisat2_index");
	make_path $cache_dir unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');

    my $done_file = "$name.done";

    my $cmd = 'nice ' . get_command_path('HISAT2_BUILD', 'hisat2-build') . " -p 32 $fasta $name && touch $done_file";

    return {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [ $fasta ],
        outputs => [
            $name . ".1.ht2",
            $name . ".2.ht2",
            $name . ".3.ht2",
            $name . ".4.ht2",
            $name . ".5.ht2",
            $name . ".6.ht2",
            $name . ".7.ht2",
            $name . ".8.ht2",
            $done_file
        ],
        description => "Indexing genome sequence with hisat2-build"
    };
}

sub create_bowtie2_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $done_files = $opts{done_files};
    my $staging_dir = $opts{staging_dir};
    my $read_type = $opts{read_type};
    my $encoding = $opts{encoding};
    my $doSeparately = $opts{doSeparately}; # for ChIP-seq pipeline, align each fastq in separate runs rather than together
    my $params = $opts{params};
    my $presets = $params->{presets} // '--sensitive';
    my $read_group = $params->{'--rg-id'} // '';

    my @tasks;

    # Add index task
    my ($index_path, $index_task) = create_bowtie_index_job($gid, $fasta);
    push @tasks, $index_task;

    # Add one or more alignment tasks
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            my $task = create_bowtie2_alignment_job(
                staging_dir => $staging_dir,
                fasta => $fasta,
                fastq => [ $file ],
                done_files => $done_files,
                index_name => $index_path,
                index_files => $index_task->{outputs},
                read_type => $read_type,
                encoding => $encoding,
                presets => $presets,
                read_group => $read_group
            );
            push @tasks, $task;
            push @sam_files, $task->{outputs}[0];
        }
    }
    else { # standard Bowtie run (all fastq's at once)
        my $task = create_bowtie2_alignment_job(
            staging_dir => $staging_dir,
            fasta => $fasta,
            fastq => $fastq,
            done_files => $done_files,
            index_name => $index_path,
            index_files => $index_task->{outputs},
            read_type => $read_type,
            encoding => $encoding,
            presets => $presets,
            read_group => $read_group
        );
        push @tasks, $task;
        push @sam_files, $task->{outputs}[0];
    }

    # Add one or more sam-to-bam tasks
    my %results;
    foreach my $file (@sam_files) {
        my $task = create_sam_to_bam_job($file, $staging_dir);
        push @tasks, $task;
        push @{$results{bam_files}}, $task->{outputs}[0];
    }

    return (\@tasks, \%results);
}

sub create_tophat_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $done_files = $opts{done_files};
    my $gff = $opts{gff};
    my $staging_dir = $opts{staging_dir};
    my $read_type = $opts{read_type};
    my $encoding = $opts{encoding};
    my $params = $opts{params};
    my $doSeparately = $opts{doSeparately}; # for ChIP-seq pipeline, align each fastq in separate runs rather than together

    my ($index, $bowtie) = create_bowtie_index_job($gid, $fasta);

    my (@tasks, @bam_files);
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            my $task = create_tophat_job(
                staging_dir => $staging_dir,
                fasta       => $fasta,
                fastq       => [ $file ],
                done_files  => $done_files,
                gff         => $gff,
                index_name  => $index,
                index_files => ($bowtie->{outputs}),
                read_type   => $read_type,
                encoding    => $encoding,
                params      => $params,
            );
            push @tasks, $task;
            push @bam_files, $task->{outputs}[0];
        }
    }
    else {
        my $task = create_tophat_job(
            staging_dir => $staging_dir,
            fasta       => $fasta,
            fastq       => $fastq,
            done_files  => $done_files,
            gff         => $gff,
            index_name  => $index,
            index_files => ($bowtie->{outputs}),
            read_type   => $read_type,
            encoding    => $encoding,
            params      => $params,
        );
        push @tasks, $task;
        push @bam_files, $task->{outputs}[0];
    }

    # Return the bam output name and jobs required
    my %results = (
        bam_files => \@bam_files
    );
    return \@tasks, \%results;
}

sub create_bowtie_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);

    my $cmd = get_command_path('BOWTIE_BUILD', 'bowtie2-build');
    my $BOWTIE_CACHE_DIR = catdir(get_genome_cache_path($gid), "bowtie_index");

    return catdir($BOWTIE_CACHE_DIR, $name), {
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fasta, 1],
            ["", $name, 0],
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($BOWTIE_CACHE_DIR, $name . ".1.bt2"),
            catfile($BOWTIE_CACHE_DIR, $name . ".2.bt2"),
            catfile($BOWTIE_CACHE_DIR, $name . ".3.bt2"),
            catfile($BOWTIE_CACHE_DIR, $name . ".4.bt2"),
            catfile($BOWTIE_CACHE_DIR, $name . ".rev.1.bt2"),
            catfile($BOWTIE_CACHE_DIR, $name . ".rev.2.bt2")
        ],
        description => "Indexing genome sequence with Bowtie"
    };
}

sub create_bowtie2_alignment_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fasta       = $opts{fasta};     # reheadered fasta file
    my $fastq       = $opts{fastq};     # array ref to fastq files
    my $done_files  = $opts{done_files}; # array ref to validation files
    my $index_name  = basename($opts{index_name});
    my $index_files = $opts{index_files};
    my $read_type   = $opts{read_type} // 'single';
    my $encoding    = $opts{encoding};
    my $presets     = $opts{presets} // '--sensitive';
    my $read_group  = $opts{read_group} // '';

    # Setup input dependencies
    my $inputs = [
        $fasta,
        @$fastq,
        @$done_files,
        @$index_files
    ];

    # Build up command/arguments string
    my $cmd = $CONF->{BOWTIE2} || 'bowtie2';
    $cmd = 'nice ' . $cmd; # run at lower priority
    $cmd .= ' -p 32';
    $cmd .= ' ' . shell_quote($presets);
    $cmd .= ' --rg-id ' . shell_quote($read_group) if $read_group;
    $cmd .= ' --phred64' if ($encoding eq '64'); # default is --phred33
    $cmd .= " -x $index_name ";

    if ($read_type eq 'paired') {
        my ($m1, $m2) = detect_paired_end(\@$fastq);
        die "error: invalid paired-end files: m1: @$m1 -- m2: @$m2" unless (@$m1 and @$m2);
        $cmd .= '-1 ' . join(',', sort @$m1) . ' -2 ' . join(',', sort @$m2);
    }
    else { # single-ended
        $cmd .= '-U ' . join(',', sort @$fastq);
    }

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.sam';
    $cmd .= " -S $output_file";

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    return {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, $output_file)
        ],
        description => "Aligning $desc using Bowtie2"
    };
}

sub create_tophat_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fasta       = $opts{fasta};
    my $fastq       = $opts{fastq};
    my $done_files  = $opts{done_files};
    my $gff         = $opts{gff};
    my $index_name  = basename($opts{index_name});
    my $index_files = $opts{index_files};
    my $read_type   = $opts{read_type} // 'single';
    my $encoding    = $opts{encoding};
    my $params      = $opts{params};

    # Optional arguments
    my $g = $params->{'-g'} // 1; #/

    # Setup input dependencies
    my $inputs = [
        $fasta,
        @$fastq,
        @$done_files,
        @$index_files
    ];
    push @$inputs, $gff if $gff;

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.bam';

    # Build up command/arguments string
    my $cmd = get_command_path('TOPHAT');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $arg_str;
    $arg_str .= $cmd . ' ';
    $arg_str .= "-G $gff " if ($gff);
    $arg_str .= "--phred64_quals " if ($encoding eq '64');
    $arg_str .= "-o . -g $g -p 32 $index_name ";

    return {
        cmd => catfile($CONF->{SCRIPTDIR}, 'tophat.pl'), # this script was created because JEX can't handle TopHat's paired-end argument syntax
        script => undef,
        args => [
            ['-read_type', $read_type, 0],
            ['-cmd_args', shell_quote($arg_str), 0],
            ['-output', $output_file, 0],
            ['-files', shell_quote(join(',', @$fastq)), 0]
        ],
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, $output_file)
        ],
        description => "Aligning sequences with TopHat"
    };
}

sub create_bismark_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $done_files = $opts{done_files};
    my $staging_dir = $opts{staging_dir};
    my $encoding = $opts{encoding};
    my $read_type = $opts{read_type};
    my $params = $opts{params};

    my ($index_path, %index_task) = create_bismark_index_job($gid, $fasta);

    my %align_task = create_bismark_alignment_job(
        staging_dir => $staging_dir,
        fasta => $fasta,
        fastq => $fastq,
        done_files => $done_files,
        index_path => $index_path,
        index_files => ($index_task{outputs}),
        encoding => $encoding,
        read_type => $read_type,
        params => $params
    );

    # Return the bam output name and jobs required
    my @tasks = ( \%index_task, \%align_task );
    my %results = (
        bam_file => $align_task{outputs}->[0]
    );
    return \@tasks, \%results;
}

sub create_bismark_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);

    my $done_file = 'bismark_genome_preparation.done';

    my $cmd = $CONF->{BISMARK_DIR} ? catfile($CONF->{BISMARK_DIR}, 'bismark_genome_preparation') : 'bismark_genome_preparation';
    my $BISMARK_CACHE_DIR = catdir(get_genome_cache_path($gid), "bismark_index");
    $cmd = "mkdir -p $BISMARK_CACHE_DIR && " .
           "cp $fasta $BISMARK_CACHE_DIR/$name.fa && " . # bismark requires fasta file to end in .fa or .fasta, not .faa
           "nice $cmd $BISMARK_CACHE_DIR && " .
           "touch $done_file";

    return $BISMARK_CACHE_DIR, (
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.1.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.2.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.3.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.4.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.rev.1.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.rev.2.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.1.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.2.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.3.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.4.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.rev.1.bt2'),
            catfile($BISMARK_CACHE_DIR, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.rev.2.bt2'),
            catfile($BISMARK_CACHE_DIR, $done_file)
        ],
        description => "Indexing genome sequence with Bismark"
    );
}

sub create_bismark_alignment_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fastq       = $opts{fastq};
    my $done_files  = $opts{done_files};
    my $index_path  = $opts{index_path};#basename($opts{index_path});
    my $index_files = $opts{index_files};
    my $encoding    = $opts{encoding};
    my $read_type   = $opts{read_type} // 'single'; #/
    my $params      = $opts{params};

    # Optional arguments
    my $N = $params->{'-N'} // 0; #/
    my $L = $params->{'-L'} // 20; #/

    # Setup input dependencies
    my $inputs = [
        @$fastq,
        @$done_files,
        @$index_files
    ];

    # Build up command/arguments string
    my $cmd = $CONF->{BISMARK_DIR} ? catfile($CONF->{BISMARK_DIR}, 'bismark') : 'bismark';
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ['-p', 4, 0], # documentation states that 4 cpus is optimal, more yields diminishing returns
        [($encoding eq '64' ? '--phred64-quals' : '--phred33-quals'), '', 0],
        ['-N', $N, 0],
        ['-L', $L, 0],
        [$index_path, '', 0]
    ];

    my ($output_bam) = @$fastq;
    $output_bam =~ s/\.fastq$//; # mdb added 7/29/16 -- remove trailing ".fastq" for cutadapt
    $output_bam =~ s/\.fq$//;    # mdb added 8/8/16  -- remove trailing ".fq" for cutadapt

    if ($read_type eq 'paired') {
        $output_bam .= '_bismark_bt2_pe.bam';
        push @$args, ['-1', shift @$fastq, 0];
        push @$args, ['-2', shift @$fastq, 0];
    }
    else { # single-ended
        $output_bam .= '_bismark_bt2.bam';
        push @$args, ['', join(' ', @$fastq), 0];
    }

    return (
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => $inputs,
        outputs => [
            $output_bam
        ],
        description => "Aligning sequences with Bismark"
    );
}

sub create_bwameth_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $done_files = $opts{done_files};
    my $staging_dir = $opts{staging_dir};
    my $read_type = $opts{read_type};
    my $params = $opts{params};

    my ($index_path, $index_task, $done_file) = create_bwameth_index_job($gid, $fasta);
    push @$done_files, $done_file;

    my %align_task = create_bwameth_alignment_job(
        staging_dir => $staging_dir,
        fasta => $fasta,
        fastq => $fastq,
        done_files => $done_files,
        index_path => $index_path,
        index_files => ($index_task->{outputs}),
        read_type => $read_type,
        params => $params
    );

    # Return the bam output name and jobs required
    my @tasks = ( $index_task, \%align_task );
    my %results = (
        bam_file => $align_task{outputs}->[0]
    );
    return \@tasks, \%results;
}

sub create_bwameth_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);

    my $done_file = 'bwameth_index.done';

    my $cmd = ($CONF->{BWAMETH} ? $CONF->{BWAMETH} : 'bwameth') . ' index';
    my $BWAMETH_CACHE_DIR = catdir(get_genome_cache_path($gid), "bwameth_index");

    $cmd = "mkdir -p $BWAMETH_CACHE_DIR && " .
           "cd $BWAMETH_CACHE_DIR && " .
           "cp $fasta . && " .
           "$cmd $name && " .
           "touch $done_file";

    return $BWAMETH_CACHE_DIR, {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($BWAMETH_CACHE_DIR, "$name.bwameth.c2t"),
            catfile($BWAMETH_CACHE_DIR, "$name.bwameth.c2t.amb"),
            catfile($BWAMETH_CACHE_DIR, "$name.bwameth.c2t.ann"),
            catfile($BWAMETH_CACHE_DIR, "$name.bwameth.c2t.bwt"),
            catfile($BWAMETH_CACHE_DIR, "$name.bwameth.c2t.pac"),
            catfile($BWAMETH_CACHE_DIR, "$name.bwameth.c2t.sa"),
            catfile($BWAMETH_CACHE_DIR, $done_file)
        ],
        description => "Indexing genome sequence with bwameth"
    },
    catfile($BWAMETH_CACHE_DIR, $done_file);
}

sub create_bwameth_alignment_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fastq       = $opts{fastq};
    my $done_files  = $opts{done_files};
    my $index_path  = $opts{index_path};#basename($opts{index_path});
    my $index_files = $opts{index_files};
    my $read_type   = $opts{read_type} // 'single'; #/
    my $params      = $opts{params};

    # Setup input dependencies
    my $inputs = [
        @$fastq,
        @$done_files,
        @$index_files
    ];

    # Build command and arguments
    my $cmd = $CONF->{BWAMETH} || 'bwameth';
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ['--reference', catfile($index_path, 'genome.faa.reheader.faa'), 0],
        ['', join(' ', @$fastq), 0],
        ['-t', 8, 0],
        ['-p', 'alignment', 0]
    ];

    return (
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, 'alignment.bam')
        ],
        description => "Aligning sequences with bwameth"
    );
}

sub create_gsnap_workflow {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $done_files = $opts{done_files};
    my $read_type  = $opts{read_type} // 'single';
    my $params = $opts{params} // {};
    my $doSeparately = $opts{doSeparately}; # for ChIP-seq pipeline, align each fastq in separate runs rather than together

    my @tasks;

    # Generate index
    my $gmap_task = create_gmap_index_job($gid, $fasta);
    push @tasks, $gmap_task;

    # Generate sam file
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            my $gsnap_task = create_gsnap_alignment_job(
                fastq => [ $file ],
                done_files => $done_files,
                gmap => $gmap_task->{outputs}->[0]->[0],
                staging_dir => $staging_dir,
                read_type => $read_type,
                params => $params
            );
            push @tasks, $gsnap_task;
            push @sam_files, $gsnap_task->{outputs}->[0];
        }
    }
    else { # standard GSNAP run (all fastq's at once)
        my $gsnap_task = create_gsnap_alignment_job(
            fastq => $fastq,
            done_files => $done_files,
            gmap => $gmap_task->{outputs}->[0]->[0],
            staging_dir => $staging_dir,
            read_type => $read_type,
            params => $params,
        );
        push @tasks, $gsnap_task;
        push @sam_files, $gsnap_task->{outputs}->[0];
    }

    # Add one or more sam-to-bam tasks
    my %results;
    foreach my $file (@sam_files) {
        # Filter sam file
        my $filter_task = create_sam_filter_job($file, $staging_dir);
        push @tasks, $filter_task;

        # Convert sam file to bam
        my $bam_task = create_sam_to_bam_job($filter_task->{outputs}->[0], $staging_dir);
        push @tasks, $bam_task;
        push @{$results{bam_files}}, $bam_task->{outputs}->[0];
    }

    return (\@tasks, \%results);
}

sub create_gmap_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);
    my $cmd = get_command_path('GMAP_BUILD');
    my $GMAP_CACHE_DIR = catdir(get_genome_cache_path($gid), "gmap_index");

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ["-D", ".", 0],
            ["-d", $name . "-index", 0],
            ["", $fasta, 1]
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            [catdir($GMAP_CACHE_DIR, $name . "-index"), 1]
        ],
        description => "Indexing genome sequence with GMAP"
    };
}

sub create_gsnap_alignment_job {
    my %opts = @_;

    # Required arguments
    my $fastq = $opts{fastq};
    my $done_files = $opts{done_files};
    my $gmap = $opts{gmap};
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $read_type = $opts{read_type} // "single";
    my $params = $opts{params};
    my $gapmode = $params->{'--gap-mode'} // "none";
    my $Q = $params->{'-Q'} // 1; #/
    my $n = $params->{'-n'} // 5; #/
    my $N = $params->{'-N'} // 1; #/
    my $nofails = $params->{'--nofails'} // 1; #/
    my $max_mismatches = $params->{'--max-mismatches'};

    my ($first_fastq) = @$fastq;
    my $output_file = basename($first_fastq) . '.sam';

    my $index_name = basename($gmap);

    my $cmd = get_command_path('GSNAP');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ["-D", ".", 0],
        ["-d", $index_name, 0],
        ["--nthreads=32", '', 0],
        ["-n", $n, 0],
        ["-N", $N, 0],
        ["--format=sam", '', 0],
        ["--gmap-mode=$gapmode", '', 1],
        ["--batch=5", '', 0]
    ];

    push @$args, ["-Q", "", 0] if $Q;
    push @$args, ["--nofails", "", 1] if $nofails;
    push @$args, ["--max-mismatches=$max_mismatches", "", 0] if $max_mismatches;
    push @$args, ['--force-single-end', '', 0] if ($read_type eq 'single');

    # Sort fastq files in case of paired-end reads,
    # see http://research-pub.gene.com/gmap/src/README
    foreach (sort @$fastq) {
        push @$args, ["", $_, 1];
    }

    push @$args, [">", $output_file, 1];

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    return {
        cmd => $cmd,
        script => undef,
# mdb removed 2/2/15 -- fails on zero-length validation input
#        options => {
#            "allow-zero-length" => JSON::false,
#        },
        args => $args,
        inputs => [
            @$fastq,
            @$done_files,
            [$gmap, 1]
        ],
        outputs => [
            catfile($staging_dir, $output_file)
        ],
        description => "Aligning $desc with GSNAP"
    };
}

sub create_bwa_workflow {
    my %opts = @_;
    my $fasta        = $opts{fasta};
    my $fastq        = $opts{fastq};
    my $done_files   = $opts{done_files};
    my $gid          = $opts{gid};
    my $encoding     = $opts{encoding};
    my $read_type    = $opts{read_type} // 'single';
    my $staging_dir  = $opts{staging_dir};
    my $params       = $opts{params} // {};
    my $doSeparately = $opts{doSeparately}; # for ChIP-seq pipeline, align each fastq in separate runs rather than together

    my @tasks;

    # Add index task
    my $index_task = create_bwa_index_job($gid, $fasta);
    push @tasks, $index_task;

    # Add one or more alignment tasks
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            my $task = create_bwa_alignment_job(
                fastq => [ $file ],
                done_files => $done_files,
                gid => $gid,
                index_files => $index_task->{outputs},
                encoding => $encoding,
                read_type => $read_type,
                staging_dir => $staging_dir,
                params => $params
            );
            push @tasks, $task;
            push @sam_files, $task->{outputs}->[0];
        }
    }
    else { # standard Bowtie run (all fastq's at once)
        my $task = create_bwa_alignment_job(
            fastq => $fastq,
            done_files => $done_files,
            gid => $gid,
            index_files => $index_task->{outputs},
            encoding => $encoding,
            read_type => $read_type,
            staging_dir => $staging_dir,
            params => $params
        );
        push @tasks, $task;
        push @sam_files, $task->{outputs}->[0];
    }

    # Add one or more sam-to-bam tasks
    my %results;
    foreach my $file (@sam_files) {
        my $task = create_sam_to_bam_job($file, $staging_dir);
        push @tasks, $task;
        push @{$results{bam_files}}, $task->{outputs}->[0];
    }

    return (\@tasks, \%results);
}

sub create_bwa_alignment_job {
    my %opts = @_;
    my $fastq       = $opts{fastq};
    my $done_files  = $opts{done_files};
    my $gid         = $opts{gid};
    my $index_files = $opts{index_files};
    #my $encoding    = $opts{encoding};
    #my $read_type   = $opts{read_type};
    my $staging_dir = $opts{staging_dir};
    my $params = $opts{params};
    my $M = $params->{'-M'} // 0;
    my $R = $params->{'-R'} // '';

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.sam';

    my $index_path = catfile(get_genome_cache_path($gid), 'bwa_index', 'genome.reheader');

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    my @args;
    push @args, ['-M', '', 0] if $M;
    push @args, ['-R', shell_quote($R), 0] if $R;
    push @args, (
        ['-t', '32',                    0],
        ['',   $index_path,             0],
        ['',   join(' ', sort @$fastq), 0],
        ['>',  $output_file,            1]
    );

	return {
        cmd => 'nice ' . get_command_path('BWA', 'bwa') . ' mem',
        args => \@args,
        inputs => [
            @$fastq,
            @$done_files,
            @$index_files
        ],
        outputs => [
            catfile($staging_dir, $output_file)
        ],
        description => "Aligning $desc using BWA-MEM"
	};
}

sub create_bwa_index_job {
    my $gid = shift;
    my $fasta = shift;

    my $cache_dir = catdir(get_genome_cache_path($gid), "bwa_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');
    my $done_file = "$name.done";

    my $cmd = 'nice ' . get_command_path('BWA', 'bwa') . " index $fasta -p $name && touch $done_file";

    return {
        cmd => $cmd,
        args => [],
        inputs => [ $fasta ],
        outputs => [
            $name . ".amb",
            $name . ".ann",
            $name . ".bwt",
            $name . ".pac",
            $name . ".sa",
            $done_file
        ],
        description => "Indexing genome sequence with BWA"
    };
}

sub generate_additional_metadata { #TODO redo arg capture in a more automated fashion
    my ($read_params, $trimming_params, $alignment_params) = @_;

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};

    if ($trimming_params && $trimming_params->{trimmer}) {
        if ($trimming_params->{trimmer} eq 'cutadapt') {
            push @annotations, 'note|cutadapt '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '-m'));
        }
        elsif ($trimming_params->{trimmer} eq 'trimgalore') {
            push @annotations, 'note|trimgalore '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--length', '-a'));
        }
    }

    switch( _get_aligner($alignment_params) ) {
        case 'hisat2'  {
            push @annotations, qq{note|hisat2_build};
            my $params = ($read_params->{encoding} eq '64' ? '--phred64' : '--phred33');
            push @annotations, 'note|hisat2 ' . $params;
        }
        case 'bowtie2' {
            my $rg = $alignment_params->{'--rg-id'};
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|bowtie2 ' . $alignment_params->{'presets'} . ($rg ? " --rg-id $rg" : '');
        }
        case 'tophat'  {
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|tophat ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-g'));
        }
        case 'bismark' {
            push @annotations, qq{note|bismark_genome_preparation};
            push @annotations, 'note|bismark ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-L'));
        }
        case 'bwameth' {
            push @annotations, qq{note|bwameth index};
            push @annotations, 'note|bwameth (default options)';
        }
        case 'bwa'     {
            my $M = $alignment_params->{'-M'};
            my $R = $alignment_params->{'-R'};
            my $args_str = ($M ? '-M' : '') . ($R ? " -R $R" : '');
            push @annotations, qq{note|bwa index};
            push @annotations, 'note|bwa mem ' . ($args_str ? $args_str : ' (default options)');
        }
        case 'gsnap'   {
            push @annotations, qq{note|gmap_build};
            push @annotations, 'note|gsnap ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-n', '-Q', '--gap-mode', '--nofails'));
        }
    }

    return \@annotations;
}

sub _get_aligner {
    my $alignment_params = shift;

    if ($alignment_params && $alignment_params->{tool}) {
        return lc($alignment_params->{tool});
    }

    return 'gsnap'; # default aligner if not specified
}

1;
