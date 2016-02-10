package CoGe::Builder::CommonTasks;

use strict;
use warnings;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use URI::Escape::JavaScript qw(escape);
use Data::Dumper;

use CoGe::Accessory::Utils qw(detect_paired_end sanitize_name to_filename to_filename_without_extension to_filename_base);
use CoGe::Accessory::IRODS qw(irods_iget irods_iput);
use CoGe::Accessory::Web qw(get_defaults split_url);
use CoGe::Core::Storage qw(get_workflow_results_file get_download_path get_sra_cache_path);
use CoGe::Core::Metadata qw(tags_to_string);

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
    send_email_job add_items_to_notebook_job create_hisat2_workflow
    export_experiment_job create_cutadapt_workflow
    create_trimgalore_job create_trimgalore_workflow
    create_bismark_alignment_job create_bismark_index_job create_bismark_workflow
    create_bwameth_alignment_job create_bwameth_index_job create_bwameth_workflow
    create_bowtie2_alignment_job create_bowtie2_workflow create_bowtie2_chipseq_workflow
);

our $CONF = CoGe::Accessory::Web::get_defaults();

sub link_results { #FIXME deprecated, remove soon
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

sub generate_results { #FIXME deprecated, remove soon
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

    return $output_file, {
        cmd  => catfile($CONF->{SCRIPTDIR}, "coge2bed.pl"),
        args => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0],
        ],
        outputs => [$output_file]
    };
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

    return $output_file, {
        cmd     => catfile($CONF->{SCRIPTDIR}, "export_NCBI_TBL.pl"),
        args    => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ["-config", $CONF->{_CONFIG_PATH}, 0]
        ],
        outputs => [$output_file]
    };
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

sub export_experiment_job {
    my %args = @_;
    my $eid = $args{eid};
    my $output = $args{output};

    return {
        cmd => catdir($CONF->{SCRIPTDIR}, "export_experiment_or_genome.pl"),
        description => "Generating experiment files",
        args => [
            ["-id", $eid, 0],
            ["-type", '"experiment"', 0],
            ["-output", $output, 1],
            ["-conf", $CONF->{_CONFIG_PATH}, 0],
            ["-dir", ".", ""]
        ],
        inputs => [],
        outputs => [$output]
    };
}

sub create_iget_job {
    my %args = @_;
    my $irods_path = $args{irods_path}; # source path
    my $local_path = $args{local_path}; # destination path

    my $dest_file = catdir($local_path, 'irods', $irods_path);
    my $done_file = $dest_file . '.done';
    my $dest_path = dirname($dest_file);
    #make_path($dest_path) unless (-r $dest_path); # mdb removed 2/9/16 -- for hypnotoad
    my $cmd;
    $cmd .= "mkdir -p $dest_path ; "; # mdb added 2/9/16 -- for hypnotoad
    $cmd .= irods_iget( $irods_path, $dest_path, { no_execute => 1 } ) . ' ; ';
    $cmd .= "touch $done_file";

    return {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [],
        outputs => [ 
            $dest_file,
            $done_file
        ],
        description => "Fetching $irods_path..."
    };
}

sub create_ftp_get_job {
    my %opts = @_;
    my $url = $opts{url};
    my $username = $opts{username} // ''; #/;
    my $password = $opts{password} // ''; #/;
    my $dest_path = $opts{dest_path};
    
    my ($filename, $path) = split_url($url);
    my $output_file = catfile($dest_path, $path, $filename);
    
    return {
        cmd => catfile($CONF->{SCRIPTDIR}, "ftp.pl"),
        script => undef,
        args => [
            ['-url',       "'".$url."'",       0],
            ['-username',  "'".$username."'",  0],
            ["-password",  "'".$password."'",  0],
            ["-dest_path", $dest_path,         0]
        ],
        inputs => [],
        outputs => [ 
            $output_file
        ],
        description => "Fetching $url..."
    };
}

sub create_data_retrieval_workflow {
    my %opts = @_;
    my $upload_dir = $opts{upload_dir};
    my $data       = $opts{data};
    
    my (@tasks, @outputs, @ncbi);
    foreach my $item (@$data) {
        my $type = lc($item->{type});
        
        # Check if the file already exists which will be the case if called
        # via the load page.
        if ($item->{path}) {
            my $filepath = catfile($upload_dir, $item->{path});
            if (-r $filepath) {
                push @outputs, $filepath;
                next;
            }
        }
        
        # Create task based on source type (IRODS, HTTP, FTP)
        my $task;
        if ($type eq 'irods') {
            my $irods_path = $item->{path};
            $irods_path =~ s/^irods//; # strip of leading "irods" from LoadExperiment page # FIXME remove this in FileSelect
            $task = create_iget_job(irods_path => $irods_path, local_path => $upload_dir);
        }
        elsif ($type eq 'http' or $type eq 'ftp') {
            $task = create_ftp_get_job(
                url => $item->{url} || $item->{path},
                username => $item->{username},
                pasword => $item->{password},
                dest_path => $upload_dir
            );
        }
        elsif ($type eq 'ncbi') {
            #TODO move file retrieval from genbank_genome_loader.pl to here
            push @ncbi, $item->{path};
        }
        elsif ($type eq 'sra') {
            my $cache_path = get_sra_cache_path();
            $task = create_fastq_dump_job(
                accn => $item->{path},
                dest_path => $cache_path
            );
        }
        
        # Add task to workflow
        if ($task) {
            push @tasks, $task;
            push @outputs, $task->{outputs}[0];
        }
    }
    
    return {
        tasks => \@tasks,
        outputs => \@outputs,
        ncbi  => \@ncbi
    };
}

sub create_fastq_dump_job {
    my %opts = @_;
    my $accn = $opts{accn};
    my $dest_path = $opts{dest_path};
    die unless $accn;
    
    my $output_file = catfile($dest_path, $accn . '.fastq');
    my $done_file = "$output_file.done";

    my $cmd = $CONF->{FASTQ_DUMP} || 'fastq-dump';

    return {
        cmd => "$cmd $accn --outdir $dest_path ; touch $done_file",
        script => undef,
        args => [],
        inputs => [],
        outputs => [
            $output_file,
            $done_file
        ],
        description => "Fetching $accn from NCBI-SRA..."
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
        add_chr => 0
    );
    @args{(keys %inputs)} = (values %inputs);

    # Check for a genome or dataset id
    return unless $args{gid};

    # Set the default basename as the id if the basename is not set
    $args{basename} = $args{gid} unless $args{basename};

    # Generate the output filename
    my $param_string = join( "-", map { $_ . $args{$_} } qw(annos cds id_type nu upa add_chr) );
    my $filename = $args{basename} . "_" . $param_string;
    if ($args{chr}) {
    	$filename .= "_" . $args{chr};
    }
    $filename .= ".gid" . $args{gid};
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
        ['-cds', $args{cds}, 0],
        ['-annos', $args{annos}, 0],
        ['-nu', $args{nu}, 0],
        ['-id_type', $args{id_type}, 0],
        ['-upa', $args{upa}, 0]
    ];
    push @$args, ['-chr', $args{chr}, 0] if (defined $args{chr});
    push @$args, ['-add_chr', $args{add_chr}, 0] if (defined $args{add_chr});
    
    # Return workflow definition
    return $output_file, {
        cmd     => catfile($CONF->{SCRIPTDIR}, "coge_gff.pl"),
        script  => undef,
        args    => $args,
        outputs => [ $output_file ],
        description => "Generating GFF..."
    };
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
    
    my @tags = ( 'VCF' ); # add VCF tag
    push @tags, @{$metadata->{tags}} if $metadata->{tags};
    my $tags_str = tags_to_string(\@tags);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', qq["$username"], 0],
            ['-name', "'".$metadata->{name}." (SNPs)". "'", 0],
            ['-desc', qq{"$desc"}, 0],
            ['-version', "'".$metadata->{version}."'", 0],
            ['-restricted', "'".$metadata->{restricted}."'", 0],
            ['-exit_without_error_for_empty_input', 1, 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', "'".$metadata->{source}."'", 0],
            ['-tags', qq{"$tags_str"}, 0],
            ['-annotations', qq["$annotations_str"], 0],
            ['-staging_dir', "./load_vcf", 0],
            ['-file_type', qq["vcf"], 0],
            ['-data_file', $vcf, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $opts->{vcf},
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            #$result_file
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
    my $name = $opts{name}; # optional name for this load
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_name = "load_experiment" . ($name ? "_$name" : '');
    my $output_path = catdir($staging_dir, $output_name);
    
    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    my $tags_str = tags_to_string($metadata->{tags});

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
            ['-annotations', qq["$annotations_str"], 0],
            ['-tags', qq["$tags_str"], 0],
            ['-staging_dir', $output_name, 0],
            ['-data_file', $input_file, 0],
            ['-normalize', $normalize, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $input_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading" . ($name ? " $name" : '') . " experiment..."
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
    
    my $cmd = 'perl ' . catfile($CONF->{SCRIPTDIR}, "load_genome.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_path = catdir($staging_dir, "load_genome");
    
    my $result_file = get_workflow_results_file($user->name, $wid);

    my $file_str = '';
    #$file_str = join(',', map { basename($_) } @$input_files) if ($input_files && @$input_files);
    $file_str = join(',', @$input_files) if ($input_files && @$input_files);
    my $irods_str = '';
    $irods_str = join(',', map { basename($_) } @$irods_files) if ($irods_files && @$irods_files);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-wid', $wid, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . ($metadata->{description} ? $metadata->{description} : '') . '"', 0],
            ['-link', '"' . ($metadata->{link} ? $metadata->{link} : '') . '"', 0],
            ['-version', '"' . $metadata->{version} . '"', 0],
            ['-restricted', ( $metadata->{restricted} ? 1 : 0 ), 0],
            ['-source_name', '"' . $metadata->{source_name} . '"', 0],
            ['-organism_id', $organism_id, 0],
            ['-type_id', '"' . ( $metadata->{type} ? $metadata->{type} : 1 ) . '"', 0], # default to "unmasked"
            ['-staging_dir', "./load_genome", 0],
            ['-fasta_files', "'".$file_str."'", 0],
            ['-irods_files', "'".$irods_str."'", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
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
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-GO', 1, 0]
    ];
    foreach (@$ncbi_accns) {
        push @$args, ['-accn', "'$_'", 0];
    }

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [],
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
            ['-config', $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
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
        ['-config', $CONF->{_CONFIG_PATH}, 0]    
    ];
    push $args, ['-nid', $nid, 0] if ($nid);

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [
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
    my $additional_metadata = $opts{additional_metadata};
    my $staging_dir = $opts{staging_dir};
    my $annotations = $opts{annotations};
    my $wid = $opts{wid};
    my $gid = $opts{gid};
    my $bam_file = $opts{bam_file};
    die unless ($user && $staging_dir && $wid && $gid && $bam_file);
    
    my $cmd = 'perl ' . catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    my $output_name = "load_bam_" . to_filename_base($bam_file);
    my $output_path = catdir($staging_dir, $output_name);
    
    my $result_file = get_workflow_results_file($user->name, $wid);
    
    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);
    
    my @tags = ( 'BAM' ); # add BAM tag
    push @tags, @{$metadata->{tags}} if $metadata->{tags};
    my $tags_str = tags_to_string(\@tags);

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', "'" . $metadata->{name} . " (BAM alignment)" . "'", 0],
            ['-desc', "'" . $metadata->{description} . "'", 0],
            ['-version', "'" . $metadata->{version} . "'", 0],
            ['-restricted', "'" . $metadata->{restricted} . "'", 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', "'" . $metadata->{source} . "'", 0],
            ['-tags', qq{"$tags_str"}, 0],
            ['-annotations', qq["$annotations_str"], 0],
            ['-staging_dir', $output_name, 0],
            ['-file_type', qq["bam"], 0],
            ['-data_file', $bam_file, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
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
    my $trimming_params = $opts{trimming_params} // {};
    my $q = $trimming_params->{'-q'} // 25;
    my $m = $trimming_params->{'-m'} // 17;
    my $read_params = $opts{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

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
    $arg_str .= "-q $q --quality-base=$encoding -m $m -o $outputs[0] ";
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

sub create_cutadapt_workflow {
    my %opts = @_;
    my $fastq1 = $opts{fastq1}; # array ref of left reads (or all reads if single-ended)
    my $fastq2 = $opts{fastq2}; # array ref of right reads (or undef if single-ended)
    my $validated = $opts{validated};
    my $staging_dir = $opts{staging_dir};
    my $read_params = $opts{read_params};
    my $trimming_params = $opts{trimming_params};
    
    my (@tasks, @outputs);

    if (not defined $fastq2) { # single-ended
        # Create cutadapt task for each file
        foreach my $file (@$fastq1) {
            my $task = create_cutadapt_job(
                fastq => $file,
                validated => $validated,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, $task->{outputs}->[0];
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
                validated => $validated,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, @{$task->{outputs}};
            push @tasks, $task;
        }
    }
    
    return ( \@tasks, \@outputs );
}

sub create_trimgalore_job {
    my %opts = @_;

    # Required params
    my $fastq = $opts{fastq};            # for single fastq file (backwards compatibility) or two paired-end fastq files (new functionality)
    my $validated = $opts{validated};    # input dependency from previous task, one or two files based on fastq arg
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
    $validated = [ $validated ] unless (ref($validated) eq 'ARRAY');
    
    my $name = join(', ', map { basename($_) } @$fastq);
    my @inputs = ( @$fastq, @$validated);

    # Build up command/arguments string
    my $cmd = $CONF->{TRIMGALORE} || 'trim_galore';
    $cmd = 'nice ' . $cmd; # run at lower priority
    
    # Create staging dir
    $cmd = "mkdir -p $staging_dir ; " . $cmd;

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

    return {
        cmd => catfile($cmd),
        script => undef,
        args => $args,
        inputs => \@inputs,
        outputs => \@outputs,
        description => "Trimming (trimgalore) $name..."
    };
}

sub create_trimgalore_workflow {
    my %opts = @_;
    my $fastq1 = $opts{fastq1} // []; # array ref of left reads (or all reads if single-ended)
    my $fastq2 = $opts{fastq2} // []; # array ref of right reads (or undef if single-ended)
    my $validated = $opts{validated};
    my $staging_dir = $opts{staging_dir};
    my $read_params = $opts{read_params};
    my $trimming_params = $opts{trimming_params};
    
    my (@tasks, @outputs);

    if ($read_params->{read_type} eq 'single') { # single-ended
        # Create trimgalore task for each file
        foreach my $file (@$fastq1, @$fastq2) {
            my $task = create_trimgalore_job(
                fastq => $file,
                validated => $validated,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, $task->{outputs}->[0];
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
                validated => $validated,
                staging_dir => $staging_dir,
                read_params => $read_params,
                trimming_params => $trimming_params
            );
            push @outputs, @{$task->{outputs}};
            push @tasks, $task;
        }
    }
    
    return ( \@tasks, \@outputs );
}

sub create_gff_generation_job {
    my %opts = @_;

    # Required params
    my $gid = $opts{gid};
    my $organism_name = $opts{organism_name};
    #my $validated = $opts{validated}; # mdb removed 1/5/15 -- why is this here?

    my $cmd = catfile($CONF->{SCRIPTDIR}, "coge_gff.pl");
    my $name = sanitize_name($organism_name) . "-1-name-0-0-id-" . $gid . "-1.gff";

    my $inputs = [];
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
            ['-config', $CONF->{_CONFIG_PATH}, 0],
        ],
        inputs => $inputs,
        outputs => [
            catdir($CONF->{CACHEDIR}, $gid, "gff", $name)
        ],
        description => "Generating genome annotations GFF file..."
    };
}

sub create_hisat2_workflow {
    my %opts = @_;
    my $fasta        = $opts{fasta};
    my $fastq        = $opts{fastq};
    my $validated    = $opts{validated};
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
                validated => $validated,
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
            validated => $validated,
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
    my $validated   = $opts{validated};
    my $gid         = $opts{gid};
    my $index_files = $opts{index_files};
    my $encoding    = $opts{encoding};
    my $read_type   = $opts{read_type};
    my $staging_dir = $opts{staging_dir};

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.sam';

	my $args = [
		['-p', '32', 0],
		['-x', catfile($CONF->{CACHEDIR}, $gid, 'hisat2_index', 'genome.reheader'), 0],
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
	
	return {
        cmd => 'nice ' . $CONF->{HISAT2},
        script => undef,
        args => $args,
        inputs => [ @$fastq, @$validated, @$index_files ],
        outputs => [ catfile($staging_dir, $output_file) ],
        description => "Aligning " . join(', ', map { to_filename_base($_) } @$fastq) . " using HISAT2..."
	};
}

sub create_hisat2_index_job {
    my $gid = shift;
    my $fasta = shift;

    my $cache_dir = catdir($CONF->{CACHEDIR}, $gid, "hisat2_index");
	make_path $cache_dir unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');
    
    return {
        cmd => 'nice ' . $CONF->{HISAT2_BUILD},
        script => undef,
        args => [
        	['-p', '32', 0],
            ['', $fasta, 0],
            ['', $name, 0],
        ],
        inputs => [ $fasta ],
        outputs => [
            $name . ".1.ht2",
            $name . ".2.ht2",
            $name . ".3.ht2",
            $name . ".4.ht2",
            $name . ".5.ht2",
            $name . ".6.ht2",
            $name . ".7.ht2",
            $name . ".8.ht2"
        ],
        description => "Indexing genome sequence with hisat2-build..."
    };
}

#sub create_bowtie2_workflow {
#    my %opts = @_;
#
#    # Required arguments
#    my $gid = $opts{gid};
#    my $fasta = $opts{fasta};
#    my $fastq = $opts{fastq};
#    my $validated = $opts{validated};
#    my $staging_dir = $opts{staging_dir};
#    my $read_type = $opts{read_type};
#    my $encoding = $opts{encoding};
#
#    my ($index_path, $index_task) = create_bowtie_index_job($gid, $fasta);
#    
#    my $align_task = create_bowtie2_alignment_job(
#        staging_dir => $staging_dir,
#        fasta => $fasta,
#        fastq => $fastq,
#        validated => $validated,
#        index_name => $index_path,
#        index_files => $index_task->{outputs},
#        read_type => $read_type,
#        encoding => $encoding
#    );
#    
#    my $bam_task = create_sam_to_bam_job($align_task->{outputs}->[0], $staging_dir);
#
#    # Return the bam output name and tasks
#    my @tasks = ( $index_task, $align_task, $bam_task );
#    my %results = ( bam_file => $bam_task->{outputs}[0] );
#    return \@tasks, \%results;
#}

sub create_bowtie2_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $validated = $opts{validated};
    my $staging_dir = $opts{staging_dir};
    my $read_type = $opts{read_type};
    my $encoding = $opts{encoding};
    my $doSeparately = $opts{doSeparately}; # for ChIP-seq pipeline, align each fastq in separate runs rather than together

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
                validated => $validated,
                index_name => $index_path,
                index_files => $index_task->{outputs},
                read_type => $read_type,
                encoding => $encoding
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
            validated => $validated,
            index_name => $index_path,
            index_files => $index_task->{outputs},
            read_type => $read_type,
            encoding => $encoding
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
    my $validated = $opts{validated};
    my $gff = $opts{gff};
    my $staging_dir = $opts{staging_dir};
    my $read_type = $opts{read_type};
    my $encoding = $opts{encoding};
    my $params = $opts{params};

    my ($index, $bowtie) = create_bowtie_index_job($gid, $fasta);
    
    my %tophat = create_tophat_job(
        staging_dir => $staging_dir,
        fasta => $fasta,
        fastq => $fastq,
        validated => $validated,
        gff => $gff,
        index_name => $index,
        index_files => ($bowtie->{outputs}),
        read_type => $read_type,
        encoding => $encoding,
        params => $params,
    );

    # Return the bam output name and jobs required
    my @tasks = ( $bowtie, \%tophat );
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
        description => "Indexing genome sequence with Bowtie..."
    };
}

sub create_bowtie2_alignment_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fasta       = $opts{fasta};     # reheadered fasta file
    my $fastq       = $opts{fastq};     # array ref to fastq files
    my $validated   = $opts{validated}; # array ref to validation files
    my $index_name  = basename($opts{index_name});
    my $index_files = $opts{index_files};
    my $read_type   = $opts{read_type} // 'single';
    my $encoding    = $opts{encoding};

    # Setup input dependencies
    my $inputs = [
        $fasta,
        @$fastq,
        @$validated,
        @$index_files
    ];

    # Build up command/arguments string
    my $cmd = $CONF->{BOWTIE2} || 'bowtie2';
    $cmd = 'nice ' . $cmd; # run at lower priority
    
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
    
    return {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, $output_file)
        ],
        description => "Aligning " . join(', ', map { to_filename_base($_) } @$fastq) . " using Bowtie2..."
    };    
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
    my $read_type   = $opts{read_type} // 'single';
    my $encoding    = $opts{encoding};
    my $params      = $opts{params};

    # Optional arguments
    my $g = $params->{'-g'} // 1; #/

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
    $arg_str .= "--phred64_quals " if ($encoding eq '64');
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
        description => "Aligning sequences with TopHat..."
    );
}

sub create_bismark_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $validated = $opts{validated};
    my $staging_dir = $opts{staging_dir};
    my $encoding = $opts{encoding};
    my $read_type = $opts{read_type};
    my $params = $opts{params};

    my ($index_path, %index_task) = create_bismark_index_job($gid, $fasta);
    
    my %align_task = create_bismark_alignment_job(
        staging_dir => $staging_dir,
        fasta => $fasta,
        fastq => $fastq,
        validated => $validated,
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
    my $BISMARK_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "bismark_index");
    $cmd = "mkdir -p $BISMARK_CACHE_DIR ; " .
           "cp $fasta $BISMARK_CACHE_DIR/$name.fa ; " . # bismark requires fasta file to end in .fa or .fasta, not .faa
           "nice $cmd $BISMARK_CACHE_DIR ; " .
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
        description => "Indexing genome sequence with Bismark..."
    );
}

sub create_bismark_alignment_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fastq       = $opts{fastq};
    my $validated   = $opts{validated};
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
        @$validated,
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
        description => "Aligning sequences with Bismark..."
    );
}

sub create_bwameth_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $validated = $opts{validated};
    my $staging_dir = $opts{staging_dir};
    my $read_type = $opts{read_type};
    my $params = $opts{params};

    my ($index_path, %index_task) = create_bwameth_index_job($gid, $fasta);
    
    my %align_task = create_bwameth_alignment_job(
        staging_dir => $staging_dir,
        fasta => $fasta,
        fastq => $fastq,
        validated => $validated,
        index_path => $index_path,
        index_files => ($index_task{outputs}),
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

sub create_bwameth_index_job {
    my $gid = shift;
    my $fasta = shift;
    my $name = to_filename($fasta);
    
    my $done_file = 'bwameth_index.done';
    
    my $cmd = ($CONF->{BWAMETH} ? $CONF->{BWAMETH} : 'bwameth') . ' index';
    my $BWAMETH_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "bwameth_index");
    
    $cmd = "mkdir -p $BWAMETH_CACHE_DIR ; " .
           "cd $BWAMETH_CACHE_DIR ; " .
           "cp $fasta . ; " . 
           "$cmd $name ; " .
           "touch $done_file";

    return $BWAMETH_CACHE_DIR, (
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
        description => "Indexing genome sequence with bwameth..."
    );
}

sub create_bwameth_alignment_job {
    my %opts = @_;

    # Required arguments
    my $staging_dir = $opts{staging_dir};
    my $fastq       = $opts{fastq};
    my $validated   = $opts{validated};
    my $index_path  = $opts{index_path};#basename($opts{index_path});
    my $index_files = $opts{index_files};
    my $read_type   = $opts{read_type} // 'single'; #/
    my $params      = $opts{params};

    # Setup input dependencies
    my $inputs = [
        @$fastq,
        @$validated,
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
        description => "Aligning sequences with bwameth..."
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
    my $bam_task = create_sam_to_bam_job($filter{outputs}->[0], $staging_dir);

    # Return the bam output name and jobs required
    my @tasks = (
        \%gmap,
        \%gsnap,
        \%filter,
        $bam_task
    );
    my %results = (
        bam_file => $bam_task->{outputs}[0]
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

sub create_sam_to_bam_job {
    my $samfile = shift;
    my $staging_dir = shift;
    
    my $filename = to_filename($samfile);
    my $cmd = $CONF->{SAMTOOLS} || 'samtools';

    return {
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
    };
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
    die unless ($input_file && $staging_dir);
    
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
    my $read_type = $params->{read_type} // "single"; #/
    my $gapmode = $params->{'--gap-mode'} // "none"; #/
    my $Q = $params->{'-Q'} // 1; #/
    my $n = $params->{'-n'} // 5; #/
    my $N = $params->{'-N'} // 1; #/
    my $nofails = $params->{'--nofails'} // 1; #/
    my $max_mismatches = $params->{'--max-mismatches'};

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
        description => "Aligning sequences with GSNAP..."
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
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-log', $log_file, 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
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
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-log', $log_file, 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
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
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-log', $log_file, 0]
    ];

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
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
    my $from = 'CoGe Support <coge.genome@gmail.com>';
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
            @$done_files
        ],
        outputs => [ 
            $done_file
        ],
        description => "Sending email..."
    };
}

1;
