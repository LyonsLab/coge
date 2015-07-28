package CoGe::Builder::CommonTasks;

use strict;
use warnings;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use URI::Escape::JavaScript qw(escape);
use Data::Dumper;

use CoGe::Accessory::Utils qw(sanitize_name to_filename);
use CoGe::Accessory::IRODS qw(irods_iget irods_iput);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_workflow_results_file get_download_path);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    generate_results link_results generate_bed 
    generate_tbl export_to_irods generate_gff generate_features copy_and_mask
    create_fasta_reheader_job create_fasta_index_job create_load_vcf_job
    create_bam_index_job create_gff_generation_job create_load_experiment_job
    create_load_genome_job create_load_genome_from_NCBI_job create_load_batch_job
    create_validate_fastq_job create_cutadapt_job create_tophat_workflow
    create_gsnap_workflow create_load_bam_job create_gunzip_job
    create_notebook_job create_bam_sort_job create_iget_job
    create_load_annotation_job create_data_retrieval_workflow
    send_email_job add_items_to_notebook_job
);

our $CONF = CoGe::Accessory::Web::get_defaults();

sub link_results {
   my ($input, $output, $result_dir, $conf) = @_;

   return {
        cmd     => catfile($CONF->{SCRIPTDIR}, "link_results.pl"),
        args    => [
            ['-input_files', escape($input), 0],
            ['-output_files', escape($output), 0],
            ['-result_dir', $result_dir, 0]
        ],
        inputs  => [$input],
        outputs => [catfile($result_dir, basename($output))],
        description => "Generating results..."
   };
}

sub generate_results {
   my ($input, $type, $result_dir, $conf, $dependency) = @_;

   return {
        cmd     => catfile($CONF->{SCRIPTDIR}, "generate_results.pl"),
        args    => [
            ['-input_files', escape($input), 0],
            ['-type', $type, 0],
            ['-result_dir', $result_dir, 0]
        ],
        inputs  => [$dependency],
        outputs => [catfile($result_dir, "1")],
        description => "Generating results..."
   };
}

sub copy_and_mask {
    my %args = @_;

    my $desc = $args{mask} ? "Copying and masking genome" : "Copying genome";
    $desc .= " (no annotations)" if $args{seq_only};
    $desc .= "...";

    my $cmd = "/copy_genome/copy_load_mask_genome.pl";

    return (
        cmd   => catfile($CONF->{SCRIPTDIR}, $cmd),
        args  => [
            ["-conf", $CONF->{_CONFIG_PATH}, 0],
            ["-gid", $args{gid}, 0],
            ["-uid", $args{uid}, 0],
            ["-wid", $args{wid}, 0],
            ["-mask", $args{mask}, 0],
            ["-staging_dir", $args{staging_dir}, 0],
            ["-sequence_only", $args{seq_only}, 0]
        ],
        description => $desc
    );
}

sub generate_bed {
    my %args = @_;

    # Check for a genome or dataset id
    return unless $args{gid};

    # Generate file name
    my $basename = $args{basename};
    my $filename = "$basename" . "_id" . $args{gid} . ".bed";
    my $path = get_download_path('genome', $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd  => catfile($CONF->{SCRIPTDIR}, "coge2bed.pl"),
        args => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0],
        ],
        outputs => [$output_file]
    );
}

sub generate_features {
    my %args = @_;

    my $filename = $args{basename} . "-gid-" . $args{gid};
    $filename .= "-prot" if $args{protein};
    $filename .= ".fasta";
    my $path = get_download_path('genome', $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd    => catfile($CONF->{SCRIPTDIR}, "export_features_by_type.pl"),
        args   => [
            ["-config", $CONF->{_CONFIG_PATH}, 0],
            ["-f", $filename, 0],
            ["-gid", $args{gid}, 0],
            ["-ftid", $args{fid}, 0],
            ["-prot", $args{protein}, 0],
        ],
        outputs => [$output_file]
    );
}

sub generate_tbl {
    my %args = @_;

    # Generate filename
    my $organism = $args{basename};
    my $filename = $organism . "id" . $args{gid} . "_tbl.txt";
    my $path = get_download_path('genome', $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd     => catfile($CONF->{SCRIPTDIR}, "export_NCBI_TBL.pl"),
        args    => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ["-config", $CONF->{_CONFIG_PATH}, 0]
        ],
        outputs => [$output_file]
    );
}

sub export_to_irods {
    my ($src, $dest, $overwrite, $done_file) = @_;

    $overwrite = 0 unless defined $overwrite;

    my $cmd = irods_iput($src, $dest, { no_execute => 1, overwrite => $overwrite });

    my $filename = basename($done_file);

   return {
        cmd => qq[$cmd && touch $filename],
        description => "Exporting file to IRODS",
        args => [],
        inputs => [$src],
        outputs => [$done_file]
   };
}

sub create_iget_job {
    my %args = @_;
    my $irods_path = $args{irods_path}; # source path
    my $local_path = $args{local_path}; # destination path

    my $dest_file = catdir($local_path, 'irods', $irods_path);
    my $dest_path = dirname($dest_file);
    make_path($dest_path) unless (-r $dest_path);
    my $cmd = irods_iget( $irods_path, $dest_path, { no_execute => 1 } );

    return {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [],
        outputs => [ $dest_file ],
        description => "Fetching $irods_path..."
    };
}

sub create_data_retrieval_workflow {
    my %opts = @_;
    my $upload_dir = $opts{upload_dir};
    my $data = $opts{data};
    
    my (@tasks, @files);
    foreach my $item (@$data) {
        my $type = lc($item->{type});
        
        # Check if the file already exists which will be the case if called
        # via the load page.  
        my $filepath = catfile($upload_dir, $item->{path});
        if (-r $filepath) {
            push @files, $filepath;
            next;
        }
        
        # Create task based on source type (IRODS, HTTP, FTP)
        my $task;
        if ($type eq 'irods') {
            my $irods_path = $item->{path};
            $irods_path =~ s/^irods//; # strip of leading "irods" from LoadExperiment page # FIXME remove this in FileSelect
            $task = create_iget_job(irods_path => $irods_path, local_path => $upload_dir);
        }
        elsif ($type eq 'http' or $type eq 'ftp') {
            #TODO
        }
        
        # Add task to workflow
        if ($task) {
            push @tasks, $task;
            push @files, $task->{outputs}[0];
        }
    }
    
    return {
        tasks => \@tasks,
        files => \@files
    };
}

sub generate_gff {
    my %inputs = @_;
    my %args = (
        annos   => 0,
        id_type => 0,
        cds     => 0,
        nu      => 0,
        upa     => 0,
    );
    @args{(keys %inputs)} = (values %inputs);

    # Check for a genome or dataset id
    return unless $args{gid};

    # Set the default basename as the id if the basename is not set
    $args{basename} = $args{gid} unless $args{basename};

    # Generate the output filename
    my @attributes = qw(annos cds id_type nu upa add_chr);
    my $param_string = join "-", map { $_ . $args{$_} } @attributes;
    my $filename = $args{basename} . "_" . $param_string;
    if ($args{chr}) {
    	$filename .= "_" . $args{chr};
    }
    $filename .= ".gid".$args{gid};
    $filename .= ".gff";
    $filename =~ s/\s+/_/g;
    $filename =~ s/\)|\(/_/g;
    my $path = get_download_path('genome', $args{gid});
    my $output_file = catfile($path, $filename);
    
    # Build argument list
    my $args = [
        ['-gid', $args{gid}, 0],
        ['-f', $filename, 0],
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        # Parameters
        ['-cds', $args{cds}, 0],
        ['-annos', $args{annos}, 0],
        ['-nu', $args{nu}, 0],
        ['-id_type', $args{id_type}, 0],
        ['-upa', $args{upa}, 0]
    ];
    push @$args, ['-chr', $args{chr}, 0] if (defined $args{chr});
    push @$args, ['-add_chr', $args{add_chr}, 0] if (defined $args{add_chr});
    
    # Return workflow definition
    return $output_file, (
        cmd     => catfile($CONF->{SCRIPTDIR}, "coge_gff.pl"),
        args    => $args,
        outputs => [$output_file],
        description => "Generating gff..."
    );
}

sub create_gunzip_job {
    my $input_file = shift;
    my $output_file = $input_file;
    $output_file =~ s/\.gz$//;

    my $cmd = $CONF->{GUNZIP} || 'gunzip';

    return {
        cmd => "$cmd -c $input_file > $output_file ;  touch $output_file.decompressed",
        script => undef,
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            "$output_file.decompressed"
        ],
        description => "Decompressing " . basename($input_file) . "..."
    };
}

sub create_fasta_reheader_job {
    my %opts = @_;

    # Required arguments
    my $fasta          = $opts{fasta};
    my $reheader_fasta = $opts{reheader_fasta} ? $opts{reheader_fasta} : to_filename($fasta) . '.reheader.faa';
    my $cache_dir      = $opts{cache_dir};

    my $cmd = catfile($CONF->{SCRIPTDIR}, "fasta_reheader.pl");

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fasta, 1],
            ["", $reheader_fasta, 0]
        ],
        inputs => [
            $fasta,
        ],
        outputs => [
            catfile($cache_dir, $reheader_fasta),
        ],
        description => "Reheader fasta file...",
    };
}

sub create_fasta_index_job {
    my %opts = @_;
    
    # Required params
    my $fasta = $opts{fasta};
    my $cache_dir = $opts{cache_dir};

    return {
        cmd => $CONF->{SAMTOOLS} || "samtools",
        script => undef,
        args => [
            ["faidx", $fasta, 1],
        ],
        inputs => [
            $fasta,
        ],
        outputs => [
            $fasta . '.fai',
        ],
        description => "Indexing fasta file...",
    };
}

sub create_bam_index_job {
    my %opts = @_;
    #print STDERR 'create_bam_index_job ', Dumper \%opts, "\n";

    # Required arguments
    my $input_file = $opts{input_file}; # bam file

    return {
        cmd => $CONF->{SAMTOOLS} || "samtools",
        script => undef,
        args => [
            ["index", $input_file, 1],
        ],
        inputs => [
            $input_file,
        ],
        outputs => [
            $input_file . '.bai'
        ],
        description => "Indexing BAM file...",
    };
}

sub create_load_vcf_job {
    my $opts = shift;

    # Required arguments
    my $metadata = $opts->{metadata};
    my $username = $opts->{username};
    my $method   = $opts->{method};
    my $staging_dir = $opts->{staging_dir};
    my $annotations = $opts->{annotations};
    my $wid = $opts->{wid};
    my $gid = $opts->{gid};
    my $vcf = $opts->{vcf};
    
    my $desc = 'Single nucleotide polymorphisms' . ($method ? " (determined by $method method)" : '');

    my $cmd = catfile(($CONF->{SCRIPTDIR}, "load_experiment.pl"));
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_vcf");
    
    my $result_file = get_workflow_results_file($username, $wid);
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', qq["$username"], 0],
            ['-name', "'".$metadata->{name}." (SNPs)". "'", 0],
            ['-desc', qq{"$desc"}, 0],
            ['-version', "'".$metadata->{version}."'", 0],
            ['-restricted', "'".$metadata->{restricted}."'", 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', "'".$metadata->{source}."'", 0],
            ['-types', qq{"SNP"}, 0],
            ['-annotations', qq["$annotations_str"], 0],
            ['-staging_dir', "./load_vcf", 0],
            ['-file_type', qq["vcf"], 0],
            ['-data_file', $vcf, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $opts->{vcf},
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading SNPs as new experiment ..."
    };
}

sub create_load_experiment_job {
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $annotations = $opts{annotations};
    my $wid = $opts{wid};
    my $gid = $opts{gid};
    my $input_file = $opts{input_file};
    my $normalize = $opts{normalize} || 0;
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_experiment");
    
    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-user_name', $user->name, 0],
            ['-name', "'" . $metadata->{name} . "'", 0],
            ['-desc', "'" . $metadata->{description} . "'", 0],
            ['-version', "'" . $metadata->{version} . "'", 0],
            ['-restricted', "'" . $metadata->{restricted} . "'", 0],
            ['-source_name', "'" . $metadata->{source_name} . "'", 0],
            #['-types', qq{"BAM"}, 0],
            ['-annotations', qq["$annotations_str"], 0],
            ['-staging_dir', "./load_experiment", 0],
            #['-file_type', qq["bam"], 0],
            ['-data_file', $input_file, 0],
            ['-normalize', $normalize, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $input_file,
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading experiment ..."
    };
}

sub create_load_genome_job {
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $wid = $opts{wid};
    my $organism_id = $opts{organism_id};
    my $input_files = $opts{input_files};
    my $irods_files = $opts{irods_files};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_genome.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_genome");
    
    my $result_file = get_workflow_results_file($user->name, $wid);

    my $file_str = join(',', map { basename($_) } @$input_files);
    my $irods_str = join(',', map { basename($_) } @$irods_files);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-wid', $wid, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . $metadata->{description} . '"', 0],
            ['-link', '"' . $metadata->{link} . '"', 0],
            ['-version', '"' . $metadata->{version} . '"', 0],
            ['-restricted', ( $metadata->{restricted} ? 1 : 0 ), 0],
            ['-source_name', '"' . $metadata->{source_name} . '"', 0],
            ['-organism_id', $organism_id, 0],
            ['-type_id', '"' . $metadata->{type} . '"', 0],
            ['-staging_dir', "./load_genome", 0],
            ['-fasta_files', "'".$file_str."'", 0],
            ['-irods_files', "'".$irods_str."'", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH}, 
            @$input_files
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading genome ..."
    };
}

sub create_load_genome_from_NCBI_job {
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $wid = $opts{wid};
    my $ncbi_accns = $opts{ncbi_accns};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "genbank_genome_loader.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_genome_from_ncbi");
    
    my $result_file = get_workflow_results_file($user->name, $wid);

    my $args = [
        ['-user_name', $user->name, 0],
        ['-wid', $wid, 0],
        ['-staging_dir', "./load_genome_from_ncbi", 0],
        ['-config', $CONF->{_CONFIG_PATH}, 1],
        ['-GO', 1, 0]
    ];
    foreach (@$ncbi_accns) {
        push @$args, ['-accn', "'$_'", 0];
    }

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [
            $CONF->{_CONFIG_PATH}
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Importing genome from NCBI ..."
    };
}

sub create_load_annotation_job {
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $wid = $opts{wid};
    my $gid = $opts{gid};
    my $input_file = $opts{input_file};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_annotation.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_annotation");
    
    my $result_file = get_workflow_results_file($user->name, $wid);
    
    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-wid', $wid, 0],
            ['-name', '"' . ($metadata->{name} ? $metadata->{name} : '') . '"', 0],
            ['-desc', '"' . ($metadata->{description} ? $metadata->{description} : ''). '"', 0],
            ['-link', '"' . ($metadata->{link} ? $metadata->{link} : '') . '"', 0],
            ['-version', '"' . $metadata->{version} . '"', 0],
            ['-source_name', '"' . $metadata->{source} . '"', 0],
            ['-gid', $gid, 0],
            ['-staging_dir', "./load_annotation", 0],
            ['-data_file', "'".$input_file."'", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH}, 
            $input_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading annotation ..."
    };
}

sub create_load_batch_job {
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $wid = $opts{wid};
    my $nid = $opts{nid};
    my $gid = $opts{gid};
    my $files = $opts{input_files};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_batch.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $file_str = join(',', @$files);

    my $args = [
        ['-user_name', $user->name, 0],
        ['-name', '"' . $metadata->{name} . '"', 0],
        ['-desc', '"' . $metadata->{description} . '"', 0],
        ['-gid', $gid, 0],
        ['-wid', $wid, 0],
        ['-staging_dir', "'".$staging_dir."'", 0],
        ['-files', "'".$file_str."'", 0],
        ['-config', $CONF->{_CONFIG_PATH}, 1]    
    ];
    push $args, ['-nid', $nid, 0] if ($nid);

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [
            $CONF->{_CONFIG_PATH}, 
            @$files
        ],
        outputs => [
            [$staging_dir, 1],
            catdir($staging_dir, 'log.done')
        ],
        description => "Loading batch experiments..."
    };
}

sub create_load_bam_job {
    my %opts = @_;
    #print STDERR "CommonTasks::create_load_bam_job ", Dumper \%opts, "\n";

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $annotations = $opts{annotations};
    my $wid = $opts{wid};
    my $gid = $opts{gid};
    my $bam_file = $opts{bam_file};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_bam");
    
    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', "'" . $metadata->{name} . " (Alignment)" . "'", 0],
            ['-desc', "'" . $metadata->{description} . "'", 0],
            ['-version', "'" . $metadata->{version} . "'", 0],
            ['-restricted', "'" . $metadata->{restricted} . "'", 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', "'" . $metadata->{source} . "'", 0],
            ['-types', qq{"BAM"}, 0],
            ['-annotations', qq["$annotations_str"], 0],
            ['-staging_dir', "./load_bam", 0],
            ['-file_type', qq["bam"], 0],
            ['-data_file', $bam_file, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $bam_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading alignment as new experiment ..."
    };
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
            ["", $fastq, 1]
        ],
        inputs => $inputs,
        outputs => [
            "$fastq.validated"
        ],
        description => "Validating " . basename($fastq) . "..."
    };
}

sub create_cutadapt_job {
    my %opts = @_;

    # Required params
    my $fastq = $opts{fastq};            # for single fastq file (backwards compatibility) or two paired-end fastq files (new functionality)
    my $validated = $opts{validated};    # input dependency from previous task, one or two files based on fastq arg
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $params = $opts{params} // {}; #/
    my $q = $params->{'-q'} // 25; #/
    my $quality = $params->{'--quality-base'} // 32; #/
    my $m = $params->{'-m'} // 17; #/
    my $read_type = $params->{read_type} // 'single'; #/

    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');
    $validated = [ $validated ] unless (ref($validated) eq 'ARRAY');
    
    my $name = join(', ', map { basename($_) } @$fastq);
    my @inputs = ( @$fastq, @$validated);
    my @outputs = map { catfile($staging_dir, to_filename($_) . '.trimmed.fastq') } @$fastq;

    # Build up command/arguments string
    my $cmd = $CONF->{CUTADAPT};
    die "ERROR: CUTADAPT is not in the config." unless $cmd;
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $arg_str;
    $arg_str .= $cmd . ' ';
    $arg_str .= "-q $q --quality-base=$quality -m $m -o $outputs[0] ";
    $arg_str .= "-p $outputs[1] " if (@$fastq > 1); # paired-end
    
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
        description => "Trimming (cutadapt) $name..."
    };
}

sub create_gff_generation_job {
    my %opts = @_;

    # Required params
    my $gid = $opts{gid};
    my $organism_name = $opts{organism_name};
    #my $validated = $opts{validated}; # mdb removed 1/5/15 -- why is this here?

    my $cmd = catfile($CONF->{SCRIPTDIR}, "coge_gff.pl");
    my $name = sanitize_name($organism_name) . "-1-name-0-0-id-" . $gid . "-1.gff";

    my $inputs = [ $CONF->{_CONFIG_PATH} ];
    #push @{$inputs}, $validated if $validated;

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-f', $name, 0],
            ['-staging_dir', '.', 0],
            ['-gid', $gid, 0],
            ['-upa', 1, 0],
            ['-id_type', "name", 0],
            ['-cds', 0, 0],
            ['-annos', 0, 0],
            ['-nu', 1, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1],
        ],
        inputs => $inputs,
        outputs => [
            catdir($CONF->{CACHEDIR}, $gid, "gff", $name)
        ],
        description => "Generating genome annotations GFF file..."
    };
}

sub create_tophat_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $validated = $opts{validated};
    my $gff = $opts{gff};
    my $staging_dir = $opts{staging_dir};
    my $params = $opts{params};

    my ($index, %bowtie) = create_bowtie_index_job($gid, $fasta);
    
    my %tophat = create_tophat_job(
        staging_dir => $staging_dir,
        fasta => $fasta,
        fastq => $fastq,
        validated => $validated,
        gff => $gff,
        index_name => $index,
        index_files => ($bowtie{outputs}),
        params => $params,
    );

    # Return the bam output name and jobs required
    my @tasks = ( \%bowtie, \%tophat );
    my %results = (
        bam_file => $tophat{outputs}->[0]
    );
    return \@tasks, \%results;
}

sub create_bowtie_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);
    
    my $cmd = $CONF->{BOWTIE_BUILD};
    my $BOWTIE_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "bowtie_index");
    die "ERROR: BOWTIE_BUILD is not in the config." unless ($cmd);

    return catdir($BOWTIE_CACHE_DIR, $name), (
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
        description => "Indexing genome sequence with Bowtie..."
    );
}

sub create_tophat_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fasta       = $opts{fasta};
    my $fastq       = $opts{fastq};
    my $validated   = $opts{validated};
    my $gff         = $opts{gff};
    my $index_name  = basename($opts{index_name});
    my $index_files = $opts{index_files};
    my $params      = $opts{params};

    # Optional arguments
    my $g = $params->{'-g'} // 1; #/
    my $read_type = $params->{read_type} // 'single'; #/

    # Setup input dependencies
    my $inputs = [
        $fasta,
        @$fastq,
        @$validated,
        @$index_files
    ];
    push @$inputs, $gff if $gff;

    # Build up command/arguments string
    my $cmd = $CONF->{TOPHAT};
    die "ERROR: TOPHAT is not in the config." unless $cmd;
    $cmd = 'nice ' . $cmd; # run at lower priority
    
    my $arg_str;
    $arg_str .= $cmd . ' ';
    $arg_str .= "-G $gff " if ($gff);
    $arg_str .= "-o . -g $g -p 32 $index_name ";

    return (
        cmd => catfile($CONF->{SCRIPTDIR}, 'tophat.pl'), # this script was created because JEX can't handle TopHat's paired-end argument syntax
        script => undef,
        args => [
            [$read_type, '', 0],
            ['"'.$arg_str.'"', '', 0],
            ['', join(' ', @$fastq), 0]
        ],
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, "accepted_hits.bam")
        ],
        description => "Aligning sequences (tophat)..."
    );
}

sub create_gsnap_workflow {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $validated = $opts{validated};
    my $params = $opts{params} // {}; #/

    # Generate index
    my %gmap = create_gmap_index_job($gid, $fasta);

    # Generate sam file
    my %gsnap = create_gsnap_job({
        fastq => $fastq,
        validated => $validated,
        gmap => $gmap{outputs}->[0]->[0],
        staging_dir => $staging_dir,
        params => $params,
    });
    
    # Filter sam file
    my %filter = create_sam_filter_job($gsnap{outputs}->[0], $staging_dir);

    # Convert sam file to bam
    my %bam = create_samtools_bam_job($filter{outputs}->[0], $staging_dir);

    # Return the bam output name and jobs required
    my @tasks = (
        \%gmap,
        \%gsnap,
        \%filter,
        \%bam
    );
    my %results = (
        bam_file => $bam{outputs}->[0]
    );
    return \@tasks, \%results;
}

sub create_gmap_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);
    my $cmd = $CONF->{GMAP_BUILD};
    my $GMAP_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "gmap_index");

    return (
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
        description => "Indexing genome sequence with GMAP..."
    );
}

sub create_samtools_bam_job {
    my $samfile = shift;
    my $staging_dir = shift;
    
    my $filename = to_filename($samfile);
    my $cmd = $CONF->{SAMTOOLS} || 'samtools';

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["view", '', 0],
            ["-bS", $samfile, 1],
            [">", $filename . ".bam", 0]
        ],
        inputs => [
            $samfile
        ],
        outputs => [
            catfile($staging_dir, $filename . ".bam")
        ],
        description => "Generating BAM file..."
    );
}

sub create_sam_filter_job {
    my ($samfile, $staging_dir) = @_;
    
    my $filename = basename($samfile);

    my $cmd = catfile($CONF->{SCRIPTDIR}, "filter_sam.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => "$cmd $filename $filename.processed",
        script => undef,
        args => [],
        inputs => [
            $samfile
        ],
        outputs => [
            catfile($staging_dir, $filename . ".processed")
        ],
        description => "Filtering SAM file..."
    );
}

sub create_bam_sort_job {
    my %opts = @_;

    # Required arguments
    my $input_file = $opts{input_file}; # bam file
    my $staging_dir = $opts{staging_dir};
    
    my $filename = to_filename($input_file);
    my $cmd = $CONF->{SAMTOOLS} || 'samtools';

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ["sort", '', 0],
            ["", $input_file, 1],
            ["", $filename . "-sorted", 1]
        ],
        inputs => [
            $input_file
        ],
        outputs => [
            catfile($staging_dir, $filename . "-sorted.bam")
        ],
        description => "Sorting BAM file..."
    };
}

sub create_gsnap_job {
    my $opts = shift;

    # Required arguments
    my $fastq = $opts->{fastq};
    my $validated = $opts->{validated};
    my $gmap = $opts->{gmap};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $params = $opts->{params};
    my $read_type = $params->{'read_type'} // "single"; #/
    my $gapmode = $params->{'--gap-mode'} // "none"; #/
    my $Q = $params->{'-Q'} // 1; #/
    my $n = $params->{'-n'} // 5; #/
    my $N = $params->{'-N'} // 1; #/
    my $nofails = $params->{'--nofails'} // 1; #/

    my $name = basename($gmap);
    
    my $cmd = $CONF->{GSNAP};
    die "ERROR: GSNAP is not in the config." unless ($cmd);
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ["-D", ".", 0],
        ["-d", $name, 0],
        ["--nthreads=32", '', 0],
        ["-n", $n, 0],
        ["-N", $N, 0],
        ["--format=sam", '', 0],
        ["--gmap-mode=$gapmode", '', 1],
        ["--batch=5", '', 0],
    ];

    push @$args, ["-Q", "", 0] if $Q;
    push @$args, ["--nofails", "", 1] if $nofails;
    push @$args, ['--force-single-end', '', 0] if ($read_type eq 'single');
    
    # Sort fastq files in case of paired-end reads, 
    # see http://research-pub.gene.com/gmap/src/README
    foreach (sort @$fastq) { 
        push @$args, ["", $_, 1];
    }
    
    push @$args, [">", $name . ".sam", 1];

    return (
        cmd => $cmd,
        script => undef,
# mdb removed 2/2/15 -- fails on zero-length validation input
#        options => {
#            "allow-zero-length" => JSON::false,
#        },
        args => $args,
        inputs => [
            @$fastq,
            @$validated,
            [$gmap, 1]
        ],
        outputs => [
            catfile($staging_dir, $name . ".sam")
        ],
        description => "Aligning sequences (gsnap)..."
    );
}

sub create_notebook_job {
    my %opts = @_;
    my $user = $opts{user};
    my $wid = $opts{wid};
    my $metadata = $opts{metadata};
    my $annotations = $opts{annotations}; # array ref
    my $staging_dir = $opts{staging_dir};
    my $done_files = $opts{done_files};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "create_notebook.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $log_file = catfile($staging_dir, "create_notebook", "log.txt");
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);
    
    my $args = [
        ['-uid', $user->id, 0],
        ['-wid', $wid, 0],
        ['-name', '"'.$metadata->{name}.'"', 0],
        ['-desc', '"'.$metadata->{description}.'"', 0],
        ['-type', 2, 0],
        ['-restricted', $metadata->{restricted}, 0],
        ['-annotations', qq{"$annotations_str"}, 0],
        ['-config', $CONF->{_CONFIG_PATH}, 1],
        ['-log', $log_file, 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
            $CONF->{_CONFIG_PATH},
            @$done_files
        ],
        outputs => [ 
            $result_file,
            $log_file
        ],
        description => "Creating notebook of results..."
    };
}


sub add_items_to_notebook_job {
    my %opts = @_;
    my $user = $opts{user};
    my $wid = $opts{wid};
    my $notebook_id = $opts{notebook_id};
    my $staging_dir = $opts{staging_dir};
    my $done_files = $opts{done_files};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "add_items_to_notebook.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $log_file = catfile($staging_dir, "add_items_to_notebook", "log.txt");
    
    my $args = [
        ['-uid', $user->id, 0],
        ['-wid', $wid, 0],
        ['-notebook_id', $notebook_id, 0],
        ['-config', $CONF->{_CONFIG_PATH}, 1],
        ['-log', $log_file, 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
            $CONF->{_CONFIG_PATH},
            @$done_files
        ],
        outputs => [ 
            $result_file,
            $log_file
        ],
        description => "Adding experiment to notebook..."
    };
}

sub add_metadata_to_results_job {
    my %opts = @_;
    my $user = $opts{user};
    my $wid = $opts{wid};
    my $annotations = $opts{annotations}; # array ref
    my $staging_dir = $opts{staging_dir};
    my $done_files = $opts{done_files};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "add_metadata_to_results.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $log_file = catfile($staging_dir, "add_metadata_to_results", "log.txt");
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);
    
    my $args = [
        ['-uid', $user->id, 0],
        ['-wid', $wid, 0],
        ['-annotations', qq{"$annotations_str"}, 0],
        ['-config', $CONF->{_CONFIG_PATH}, 1],
        ['-log', $log_file, 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
            $CONF->{_CONFIG_PATH},
            @$done_files
        ],
        outputs => [ 
            $result_file,
            $log_file
        ],
        description => "Adding metadata to results..."
    };
}

sub send_email_job {
    my %opts = @_;
    my $from = $opts{from};
    my $to = $opts{to};
    my $subject = $opts{subject};
    my $body = $opts{body};
    my $done_files = $opts{done_files};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "send_email.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $staging_dir = $opts{staging_dir};
    my $done_file = catfile($staging_dir, "send_email.done");
    
    my $args = [
        ['-from', '"'.escape($from).'"', 0],
        ['-to', '"'.escape($to).'"', 0],
        ['-subject', '"'.escape($subject).'"', 0],
        ['-body', '"'.escape($body).'"', 0],
        ['-done_file', '"'.$done_file.'"', 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
            $CONF->{_CONFIG_PATH},
            @$done_files
        ],
        outputs => [ 
            $done_file
        ],
        description => "Sending email..."
    };
}

1;
