package CoGe::Builder::CommonTasks;

use strict;
use warnings;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use JSON qw(encode_json);
use URI::Escape::JavaScript qw(escape);
use Data::Dumper;
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils qw(detect_paired_end sanitize_name to_filename to_filename_without_extension to_filename_base);
use CoGe::Accessory::IRODS qw(irods_iget irods_iput irods_set_env);
use CoGe::Accessory::Web qw(get_defaults get_command_path split_url);
use CoGe::Core::Storage qw(get_workflow_results_file get_download_path get_sra_cache_path get_genome_cache_path);
use CoGe::Core::Metadata qw(tags_to_string);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    export_to_irods generate_gff generate_features copy_and_mask
    create_fasta_reheader_job create_fasta_index_job create_load_vcf_job create_sam_to_bam_job
    create_bam_index_job create_gff_generation_job create_load_experiment_job create_load_batch_job
    create_load_bam_job create_gunzip_job create_notebook_job create_bam_sort_job create_iget_job
    create_load_annotation_job create_data_retrieval_workflow create_sam_filter_job
    export_experiment_job create_bgzip_job create_tabix_index_job create_sumstats_job
    add_workflow_result create_image_job add_metadata_to_results_job create_bigwig_to_wig_job
    create_transdecoder_longorfs_job create_transdecoder_predict_job
);

our $CONF = CoGe::Accessory::Web::get_defaults();

# mdb removed 10/7/16 -- deprecated
#sub link_results {
#   my ($input, $output, $result_dir, $conf) = @_;
#
#   return {
#        cmd     => catfile($CONF->{SCRIPTDIR}, "link_results.pl"),
#        args    => [
#            ['-input_files', escape($input), 0],
#            ['-output_files', escape($output), 0],
#            ['-result_dir', $result_dir, 0]
#        ],
#        inputs  => [$input],
#        outputs => [catfile($result_dir, basename($output))],
#        description => "Generating results..."
#   };
#}

# mdb removed 10/7/16 -- deprecated
#sub generate_results {
#   my ($input, $type, $result_dir, $conf, $dependency) = @_;
#
#   return {
#        cmd     => catfile($CONF->{SCRIPTDIR}, "generate_results.pl"),
#        args    => [
#            ['-input_files', escape($input), 0],
#            ['-type', $type, 0],
#            ['-result_dir', $result_dir, 0]
#        ],
#        inputs  => [$dependency],
#        outputs => [catfile($result_dir, "1")],
#        description => "Generating results..."
#   };
#}

sub add_workflow_result {
    my %opts = @_;
    my $username = $opts{username};
    my $wid = $opts{wid};
    my $result = encode_json($opts{result});
    my $dependency = $opts{dependency};
    
    my $result_file = get_workflow_results_file($username, $wid);

    return {
        cmd  => catfile($CONF->{SCRIPTDIR}, "add_workflow_result.pl"),
        args => [
            ['-user_name', $username, 0],
            ['-wid', $wid, 0],
            ['-result', "'".$result."'", 0]
        ],
        inputs  => [$dependency],
        outputs => [$result_file],
        description => "Adding workflow result"
    };
}

sub copy_and_mask {
    my %args = @_;

    my $desc = $args{mask} ? "Copying and masking genome" : "Copying genome";
    $desc .= " (no annotations)" if $args{seq_only};

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

sub generate_features {
    my %args = @_;

    my $filename = $args{basename} . "-gid-" . $args{gid};
    $filename .= "-prot" if $args{protein};
    $filename .= ".fasta";
    my $path = get_genome_cache_path($args{gid});
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

sub export_to_irods {
    my ($src, $dest, $overwrite, $done_file) = @_;

    irods_set_env(catfile($CONF->{_HOME_PATH}, 'irodsEnv')); # mdb added 2/9/16 -- for hypnotoad, use www-data's irodsEnvFile
    my $cmd = irods_iput($src, $dest, { no_execute => 1, overwrite => ($overwrite // 0) });

    my $filename = basename($done_file);

    return {
        cmd => qq[$cmd && touch $done_file],
        description => "Exporting file to IRODS " . $dest,
        args => [],
        inputs => [ $src ],
        outputs => [ $done_file ]
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
    $cmd .= "mkdir -p $dest_path && "; # mdb added 2/9/16 -- for hypnotoad
    irods_set_env(catfile($CONF->{_HOME_PATH}, 'irodsEnv')); # mdb added 2/9/16 -- for hypnotoad, use www-data's irodsEnvFile
    $cmd .= irods_iget( $irods_path, $dest_path, { no_execute => 1 } ) . ' && ';
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
        description => "Fetching $irods_path"
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
        description => "Fetching $url"
    };
}

sub create_data_retrieval_workflow {
    my %opts = @_;
    my $upload_dir = $opts{upload_dir};
    my $data       = $opts{data};
    
    my (@tasks, @outputs, @ncbi);
    foreach my $item (@$data) {
        my $type = lc($item->{type} // '');
        
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
        cmd => "mkdir -p $dest_path && $cmd $accn --outdir $dest_path && touch $done_file",
        script => undef,
        args => [],
        inputs => [],
        outputs => [
            $output_file,
            $done_file
        ],
        description => "Fetching $accn from NCBI-SRA"
    };
}

#FIXME dup'ed below -- mdb 11/7/16
sub create_gff_generation_job {
    my %opts = @_;

    # Required params
    my $gid = $opts{gid};
    my $organism_name = $opts{organism_name};
    #my $validated = $opts{validated}; # mdb removed 1/5/15 -- why is this here?

    my $cmd = catfile($CONF->{SCRIPTDIR}, "coge_gff.pl");
    my $name = sanitize_name($organism_name) . "-1-name-0-0-id-" . $gid . "-1.gff";

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-f', $name, 1],
            ['-staging_dir', '.', 0],
            ['-gid', $gid, 0],
            ['-upa', 1, 0],
            ['-id_type', "name", 0],
            ['-cds', 0, 0],
            ['-annos', 0, 0],
            ['-nu', 1, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0],
        ],
        inputs => [],
        outputs => [
            catfile(get_genome_cache_path($gid), $name)
        ],
        description => "Generating genome annotations GFF file"
    };
}

#FIXME dup'ed above -- mdb 11/7/16
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
    my $path = get_genome_cache_path($args{gid});
    my $output_file = catfile($path, $filename);
    
    # Build argument list
    my $args = [
        ['-gid', $args{gid}, 0],
        ['-f', $filename, 1],
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-cds', $args{cds}, 0],
        ['-annos', $args{annos}, 0],
        ['-nu', $args{nu}, 0],
        ['-id_type', $args{id_type}, 0],
        ['-upa', $args{upa}, 0]
    ];
    push @$args, ['-chr', $args{chr}, 0] if (defined $args{chr});
    push @$args, ['-add_chr', $args{add_chr}, 0] if (defined $args{add_chr});
    
    return $output_file, {
        cmd     => catfile($CONF->{SCRIPTDIR}, "coge_gff.pl"),
        script  => undef,
        args    => $args,
        inputs  => [],
        outputs => [ $output_file ],
        description => "Generating GFF"
    };
}

sub create_image_job {
    my $input_file = shift;
    my $staging_dir = shift;
    
    return {
        cmd => catfile($CONF->{SCRIPTDIR}, 'create_image.pl'),
        script => undef,
        args => [
            ['', $input_file, 0],
            ['', $CONF->{_CONFIG_PATH}, 0],
            ['', "$input_file.log", 0]
        ],
        inputs => [
            $input_file
        ],
        outputs => [
            catfile($staging_dir, "$input_file.log")
        ],
        description => "Loading image into database"
    };
}

sub create_gunzip_job {
    my $input_file = shift;
    my $output_file = $input_file;
    $output_file =~ s/\.gz$//;

    my $cmd = get_command_path('GUNZIP');

    return {
        cmd => "$cmd -c $input_file > $output_file && touch $output_file.decompressed",
        script => undef,
        args => [],
        inputs => [
            $input_file,
            $input_file . '.done' # ensure file is done transferring
        ],
        outputs => [
            $output_file,
            "$output_file.decompressed"
        ],
        description => "Decompressing " . basename($input_file)
    };
}

sub create_bgzip_job {
    my $input_file = shift;
    my $output_file = $input_file . '.bgz';

    my $cmd = get_command_path('BGZIP');

    return {
        cmd => "$cmd -c $input_file > $output_file && touch $output_file.done",
        script => undef,
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            "$output_file.done"
        ],
        description => "Compressing " . basename($input_file) . " with bgzip"
    };
}

sub create_tabix_index_job {
    my $input_file = shift;
    my $index_type = shift; 
    my $output_file = $input_file . '.tbi';

    my $cmd = $CONF->{TABIX} || 'tabix';

    return {
        cmd => "$cmd -p $index_type $input_file && touch $output_file.done",
        script => undef,
        args => [],
        inputs => [
            $input_file,
            $input_file . '.done'
        ],
        outputs => [
            $output_file,
            "$output_file.done"
        ],
        description => "Indexing " . basename($input_file)
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
        description => "Reheader fasta file",
    };
}

sub create_fasta_index_job {
    my %opts = @_;
    
    # Required params
    my $fasta = $opts{fasta};

    return {
        cmd => get_command_path('SAMTOOLS'),
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
        description => "Indexing FASTA file",
    };
}

sub create_bam_index_job {
    my %opts = @_;

    # Required arguments
    my $input_file = $opts{input_file}; # bam file

    return {
        cmd => get_command_path('SAMTOOLS'),
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
        description => "Indexing BAM file",
    };
}

sub create_bigwig_to_wig_job {
    my %opts = @_;
    my $input_file  = $opts{input_file};
    my $staging_dir = $opts{staging_dir};
    
    my $filename = basename($input_file) . '.wig';
    my $output_file = catfile($staging_dir, $filename);
    my $done_file = $output_file . '.done';
    
    my $cmd = $CONF->{BIGWIGTOWIG} || 'bigWigToWig';

    return {
        cmd => "mkdir -p $staging_dir && $cmd $input_file $output_file && touch $done_file",
        script => undef,
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            $done_file
        ],
        description => 'Converting BigWig to WIG format'
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

    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
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
            ['-user_name',   qq["$username"], 0],
            ['-name',        ($metadata->{name} ? shell_quote($metadata->{name}." (SNPs)") : '""'), 0],
            ['-desc',        ($desc ? shell_quote($desc) : '""'), 0],
            ['-version',     ($metadata->{version} ? shell_quote($metadata->{version}) : '""'), 0],
            ['-restricted',  shell_quote($metadata->{restricted}), 0],
            ['-exit_without_error_for_empty_input', 1, 0],
            ['-gid',         $gid, 0],
            ['-wid',         $wid, 0],
            ['-source_name', ($metadata->{source_name} ? shell_quote($metadata->{source_name}) : '""'), 0],
            ['-tags',        shell_quote($tags_str), 0],
            ['-annotations', shell_quote($annotations_str), 0],
            ['-staging_dir', "./load_vcf", 0],
            ['-file_type',   "vcf", 0],
            ['-data_file',   $vcf, 0],
            ['-config',      $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $opts->{vcf},
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            #$result_file
        ],
        description => "Loading SNPs as new experiment"
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
            ['-gid',         $gid, 0],
            ['-wid',         $wid, 0],
            ['-user_name',   $user->name, 0],
            ['-name',        ($metadata->{name} ? shell_quote($metadata->{name}) : '""'), 0],
            ['-desc',        ($metadata->{description} ? shell_quote($metadata->{description}) : '""'), 0],
            ['-version',     ($metadata->{version} ? shell_quote($metadata->{version}) : '""'), 0],
            ['-restricted',  shell_quote($metadata->{restricted}), 0],
            ['-source_name', ($metadata->{source_name} ? shell_quote($metadata->{source_name}) : '""'), 0],
            ['-annotations', ($annotations_str ? shell_quote($annotations_str) : '""'), 0],
            ['-tags',        shell_quote($tags_str), 0],
            ['-staging_dir', $output_name, 0],
            ['-data_file',   $input_file, 0],
            ['-normalize',   $normalize, 0],
            ['-disable_range_check', '', 0], # mdb added 8/26/16 COGE-270 - allow values outside of [-1, 1]
            ['-config',      $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $input_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading" . ($name ? " $name" : '') . " experiment"
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
            ['-name', ($metadata->{name} ? shell_quote($metadata->{name}) : '""'), 0],
            ['-desc', ($metadata->{description} ? shell_quote($metadata->{description}) : '""'), 0],
            ['-link', ($metadata->{link} ? shell_quote($metadata->{link}) : '""'), 0],
            ['-version', ($metadata->{version} ? shell_quote($metadata->{version}) : '""'), 0],
            ['-source_name', ($metadata->{source} ? shell_quote($metadata->{source}) : '""'), 0],
            ['-gid', $gid, 0],
            ['-staging_dir', "./load_annotation", 0],
            ['-data_file', shell_quote($input_file), 0],
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
        description => "Loading annotation"
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
        ['-name', shell_quote($metadata->{name}), 0],
        ['-desc', shell_quote($metadata->{description}), 0],
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
        description => "Loading batch experiments"
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
            ['-user_name',   $user->name, 0],
            ['-name',        ($metadata->{name} ? shell_quote($metadata->{name} . " (BAM alignment)") : '""'), 0],
            ['-desc',        shell_quote($metadata->{description}), 0],
            ['-version',     shell_quote($metadata->{version}), 0],
            ['-restricted',  shell_quote($metadata->{restricted}), 0],
            ['-gid',         $gid, 0],
            ['-wid',         $wid, 0],
            ['-source_name', shell_quote($metadata->{source_name}), 0],
            ['-tags',        shell_quote($tags_str), 0],
            ['-annotations', shell_quote($annotations_str), 0],
            ['-staging_dir', $output_name, 0],
            ['-file_type',   'bam', 0],
            ['-data_file',   $bam_file, 0],
            ['-config',      $CONF->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $bam_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading alignment as new experiment"
    };
}

sub create_sam_to_bam_job {
    my $samfile = shift;
    my $staging_dir = shift;
    
    my $filename = to_filename($samfile);
    my $cmd = get_command_path('SAMTOOLS');

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
        description => "Converting SAM file to BAM"
    };
}

sub create_sam_filter_job {
    my ($samfile, $staging_dir) = @_;
    
    my $filename = basename($samfile);

    my $cmd = catfile($CONF->{SCRIPTDIR}, "filter_sam.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return {
        cmd => "$cmd $filename $filename.processed",
        script => undef,
        args => [],
        inputs => [
            $samfile
        ],
        outputs => [
            catfile($staging_dir, $filename . ".processed")
        ],
        description => "Filtering SAM file"
    };
}

sub create_bam_sort_job {
    my %opts = @_;

    # Required arguments
    my $input_file = $opts{input_file}; # bam file
    my $staging_dir = $opts{staging_dir};
    die unless ($input_file && $staging_dir);
    
    my $filename = to_filename($input_file);
    my $cmd = get_command_path('SAMTOOLS');

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
        description => "Sorting BAM file"
    };
}

sub create_sumstats_job {
    my %opts = @_;

    # Required arguments
    my $vcf = $opts{vcf};
    my $gff = $opts{gff};
    my $fasta = $opts{fasta};
    my $output_path = $opts{output_path};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "popgen/sumstats.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;
    
    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-vcf',    $vcf,         0],
            ['-gff',    $gff,         0],
            ['-fasta',  $fasta,       0],
            ['-output', $output_path, 0],
            ['-debug',  '',           0]
        ],
        inputs => [
            $vcf,
            $gff,
            $fasta
        ],
        outputs => [
            catfile($output_path, "sumstats.done"),
        ],
        description => "Calculating summary statistics"
    };
}

sub create_transdecoder_longorfs_job {
    my $input_file = shift;

    my $cmd = catfile($CONF->{TRANSDECODER}, 'TransDecoder.LongOrfs');
    
    my $done_file = $input_file . '.longorfs.done';
    my $output_path = $input_file . '.transdecoder_dir';

    return {
        cmd => "$cmd -t $input_file && touch $done_file",
        script => undef,
        args => [],
        inputs => [
            $input_file,
        ],
        outputs => [
            [$output_path, '1'],
            $done_file
        ],
        description => "Running TransDecoder.LongOrfs"
    };
}

sub create_transdecoder_predict_job {
    my $input_file = shift;
    my $input_dir = shift; 
    my $dependency_file = shift;

    my $cmd = catfile($CONF->{TRANSDECODER}, 'TransDecoder.Predict');
    
    my $done_file = $input_file . '.predict.done';
    my $output_file = $input_file . '.transdecoder.gff3';

    return {
        cmd => "cd $input_dir && $cmd -t $input_file && touch $done_file", # transdecoder won't work unless run from the dir that contains the input dir
        script => undef,
        args => [],
        inputs => [
            $input_file,
            $dependency_file
        ],
        outputs => [
            $output_file,
            $done_file
        ],
        description => "Running TransDecoder.Predict"
    };
}

sub add_metadata_to_results_job {
    my %opts = @_;
    my $user = $opts{user};
    my $wid = $opts{wid};
    my $item_id = $opts{item_id};
    my $item_type = $opts{item_type};
    my $annotations = $opts{annotations}; # array ref
    my $locked = $opts{locked} // 1;
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
        ['-locked', $locked, 0],
        ['-annotations', qq{"$annotations_str"}, 0],
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-log', $log_file, 0]
    ];
    
    if ($item_id && $item_type) {
        push @$args, ['-item_id', $item_id, 0];
        push @$args, ['-item_type', $item_type, 0];  
    }

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => [ 
            @$done_files
        ],
        outputs => [ 
            $log_file
        ],
        description => "Adding metadata to results"
    };
}

1;
