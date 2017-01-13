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
use CoGe::Accessory::TDS;
use CoGe::Core::Storage qw(get_workflow_results_file get_download_path get_sra_cache_path get_genome_cache_path);
use CoGe::Core::Metadata qw(tags_to_string);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    export_to_irods generate_gff generate_features copy_and_mask
    create_fasta_reheader_job create_fasta_index_job create_load_vcf_job
    create_bam_index_job create_gff_generation_job create_load_experiment_job create_load_batch_job
    create_load_bam_job create_notebook_job create_bam_sort_job create_iget_job
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

sub create_data_retrieval_workflow { #TODO remove ASAP: replace with DataRetrieval.pm
    my %opts = @_;
    my $upload_dir = $opts{upload_dir};
    my $data       = $opts{data};
    my $params     = $opts{params};
    
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
                dest_path => $cache_path,
                read_params => $params->{read_params}
            );
            if ($task) {
                push @tasks, $task;
                my @fastq_outputs = grep { $_ =~ /\.fastq$/ } @{$task->{outputs}}; # kludge
                push @outputs, @fastq_outputs;
                next;
            }
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
    my $read_params = $opts{read_params};
    my $read_type = $read_params->{read_type} // 'single';

    my $cmd = $CONF->{FASTQ_DUMP} || 'fastq-dump';
    $cmd .= ' --split-files' if ($read_type eq 'paired');

    my $output_filepath = catfile($dest_path, $accn);

    my (@output_files, @done_files);
    if ($read_type eq 'paired') {
        @output_files = (
            $output_filepath . '_1.fastq',
            $output_filepath . '_2.fastq'
        );
        @done_files = (
            $output_filepath . '_1.fastq.done',
            $output_filepath . '_2.fastq.done'
        );
    }
    else {
        @output_files = ( $output_filepath . '.fastq');
        @done_files   = ( $output_filepath . '.fastq.done' );
    }

    return {
        cmd => "mkdir -p $dest_path && $cmd $accn --outdir $dest_path && touch " . join(' ', @done_files),
        script => undef,
        args => [],
        inputs => [],
        outputs => [
            @output_files,
            @done_files
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
