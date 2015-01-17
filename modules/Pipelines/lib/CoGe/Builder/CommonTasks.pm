package CoGe::Builder::CommonTasks;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);
use URI::Escape::JavaScript qw(escape);
use Data::Dumper;

use CoGe::Accessory::Utils qw(sanitize_name to_filename);
use CoGe::Accessory::IRODS qw(irods_iput);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Genome qw(get_download_path);
use CoGe::Core::Storage qw(get_genome_file);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    generate_results link_results generate_bed export_experiment
    generate_tbl export_to_irods generate_gff generate_features copy_and_mask
    create_fasta_reheader_job create_fasta_index_job create_load_vcf_job
    create_bam_index_job create_gff_generation_job
    create_alignment_workflow create_load_experiment_job
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
    my %args = (
        mask => 0,
        seq_only => 0,
        @_
    );

    my $desc = $args{mask} ? "Copying and masking genome" : "Copying genome";
    $desc .= " (no annotations)" if $args{seq_only};
    $desc .= "...";

    my $cmd = "/copy_genome/copy_load_mask_genome.pl";

    return (
        cmd   => catfile($args{script_dir}, $cmd),
        args  => [
            ["-conf", $args{conf}, 0],
            ["-gid", $args{gid}, 0],
            ["-uid", $args{uid}, 0],
            ["-mask", $args{mask}, 0],
            ["-staging_dir", $args{staging_dir}, 0],
            ["-result_dir", $args{result_dir}, 0],
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
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd  => catfile($args{script_dir}, "coge2bed.pl"),
        args => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $args{conf}, 0],
        ],
        outputs => [$output_file]
    );
}

sub generate_features {
    my %args = @_;

    my $filename = $args{basename} . "-gid-" . $args{gid};
    $filename .= "-prot" if $args{protein};
    $filename .= ".fasta";
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd    => catfile($args{script_dir}, "export_features_by_type.pl"),
        args   => [
            ["-config", $args{conf}, 0],
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
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd     => catfile($args{script_dir}, "export_NCBI_TBL.pl"),
        args    => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ["-config", $args{conf}, 0]
        ],
        outputs => [$output_file]
    );
}

sub export_experiment {
    my ($params, $output, $conf) = @_;

    return (
        cmd => catdir($CONF->{SCRIPTDIR}, "export_experiment.pl"),
        description => "Generating experiment files",
        args => [
            ["-eid", $params->{eid}, 0],
            ["-output", $output, 1],
            ["-conf", $CONF->{_CONFIG_PATH}, 0],
            ["-dir", ".", ""]
        ],
        inputs => [],
        outputs => [$output]
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

sub generate_gff {
    my ($inputs, $conf) = @_;

    my %args = (
        annos   => 0,
        id_type => 0,
        cds     => 0,
        nu      => 0,
        upa     => 0,
    );

    @args{(keys $inputs)} = (values $inputs);

    # Check for a genome or dataset id
    return unless $args{gid};

    # Set the default basename as the id if the basename is not set
    $args{basename} = $args{gid} unless $args{basename};

    # Generate the output filename
    my $organism = "gff";
    my @attributes = qw(annos cds id_type nu upa);
    my $param_string = join "-", map { $_ . $args{$_} } @attributes;
    my $filename = $args{basename} . "_" . $param_string . ".gff";
    $filename =~ s/\s+/_/g;
    $filename =~ s/\)|\(/_/g;
    my $path = get_download_path($CONF->{SECTEMPDIR}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd     => catfile($CONF->{SCRIPTDIR}, "coge_gff.pl"),
        args    => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0],
            # Parameters
            ['-cds', $args{cds}, 0],
            ['-annos', $args{annos}, 0],
            ['-nu', $args{nu}, 0],
            ['-id_type', $args{id_type}, 0],
            ['-upa', $args{upa}, 0],
        ],
        outputs => [$output_file],
        description => "Generating gff..."
    );
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
        description => "Index fasta file...",
    };
}

sub create_bam_index_job {
    my %opts = @_;

    # Required arguments
    my $input_bam = $opts{input_bam};

    return {
        cmd => $CONF->{SAMTOOLS} || "samtools",
        script => undef,
        args => [
            ["index", $input_bam, 1],
        ],
        inputs => [
            $input_bam,
        ],
        outputs => [
            qq[$input_bam.bai]
        ],
        description => "Index bam file...",
    };
}

sub create_load_vcf_job {
    my $opts = shift;

    # Required arguments
    my $metadata = $opts{metadata};
    my $username = $opts->{username};
    my $staging_dir = $opts->{staging_dir};
    my $result_dir = $opts->{result_dir};
    my $annotations = $opts->{annotations};
    my $wid = $opts->{wid};
    my $gid = $opts->{gid};
    my $vcf = $opts->{vcf};

    my $cmd = catfile(($CONF->{SCRIPTDIR}, "load_experiment.pl"));
    my $output_path = catdir($staging_dir, "load_vcf");

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', qq["$username"], 0],
            ['-name', "'".$metadata->{name}." (SNPs)". "'", 0],
            ['-desc', qq{"Single nucleotide polymorphisms"}, 0],
            ['-version', "'".$metadata->{version}."'", 0],
            ['-restricted', "'".$metadata->{restricted}."'", 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', "'".$metadata->{source}."'", 0],
            ['-types', qq{"SNP"}, 0],
            ['-annotations', qq["$annotations"], 0],
            ['-staging_dir', "./load_vcf", 0],
            ['-file_type', qq["vcf"], 0],
            ['-result_file', catfile($result_dir, 'snps'), 0],
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
        ],
        description => "Load SNPs as new experiment ..."
    };
}

sub create_load_experiment_job {
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $result_dir = $opts{result_dir};
    my $annotations = $opts{annotations} || '';
    my $wid = $opts{wid};
    my $gid = $opts{gid};
    my $input_file = $opts{input_file};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    my $output_path = catdir($staging_dir, "load_experiment");

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
            ['-source_name', "'" . $metadata->{source} . "'", 0],
            #['-types', qq{"BAM"}, 0],
            ['-annotations', qq["$annotations"], 0],
            ['-staging_dir', "./load_experiment", 0],
            #['-file_type', qq["bam"], 0],
            ['-result_file', catfile($result_dir, 'experiment'), 0],
            ['-data_file', $input_file, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $input_file,
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
        ],
        description => "Loading experiment ..."
    };
}

sub create_load_bam_job {
    my %opts = @_;
    #print STDERR "CommonTasks::create_load_bam_job ", Dumper \%opts, "\n";

    # Required arguments
    my $user = $opts{user};
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $result_dir = $opts{result_dir};
    my $annotations = $opts{annotations};
    my $wid = $opts{wid};
    my $gid = $opts{gid};
    my $bam_file = $opts{bam_file};
    
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    my $output_path = catdir($staging_dir, "load_bam");

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
            ['-annotations', qq["$annotations"], 0],
            ['-staging_dir', "./load_bam", 0],
            ['-file_type', qq["bam"], 0],
            ['-result_file', catfile($result_dir, 'mapped_reads'), 0],
            ['-data_file', $bam_file, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $bam_file,
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
        ],
        description => "Loading alignment as new experiment ..."
    };
}

sub create_validate_fastq_job {
    my $fastq = shift;

    my $cmd = catfile($CONF->{SCRIPTDIR}, "validate_fastq.pl");

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fastq, 1]
        ],
        inputs => [
            $fastq
        ],
        outputs => [
            "$fastq.validated"
        ],
        description => "Validating fastq file..."
    };
}

sub create_cutadapt_job {
    my %opts = @_;

    # Required params
    my $fastq = $opts{fastq};
    my $validated = $opts{validated};
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $cutadapt_params = $opts->{cutadapt_params} // {}; #/
    my $q = $cutadapt_params->{q} // 25; #/
    my $quality = $cutadapt_params->{quality} // 64; #/
    my $m = $cutadapt_params->{m} // 17; #/

    my $inputs = [ $fastq ];
    push @{$inputs}, $validated if $validated;

    my $name = to_filename($fastq);
    my $cmd = $CONF->{CUTADAPT};

    return {
        cmd => qq[$cmd > /dev/null],
        script => undef,
        args => [
            ['-q', $q, 0],
            ["--quality-base=$quality", '', 0],
            ["-m", $m, 0],
            ['', $fastq, 1],
            ['-o', $name . '.trimmed.fastq', 1],
        ],
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, $name . '.trimmed.fastq')
        ],
        description => "Running cutadapt..."
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
        description => "Generating genome annotations gff..."
    };
}

sub create_tophat_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $gff = $opts{gff};
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $alignment_params = $opts{alignment_params} // {}; #/

    my ($index, %bowtie) = create_bowtie_index_job($gid, $fasta);
    my %tophat = create_tophat_job({
        staging_dir => $staging_dir,
        fastq => $fastq,
        fasta => $fasta,
        gff   => $gff,
        index_name => $index,
        index_files => ($bowtie{outputs}),
        g => $alignment_params->{g},
    });

    # Return the bam output name and jobs required
    return @{$tophat{outputs}}[0], (
        \%bowtie,
        \%tophat
    );
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
    my $opts = shift;

    # Required arguments
    my $fasta = $opts->{fasta};
    my $fastq = $opts->{fastq};
    my $gff = $opts->{gff};
    my $staging_dir = $opts->{staging_dir};
    my @index_files = @{$opts->{index_files}};
    my $name = basename($opts->{index_name});

    # Optional arguments
    my $g = $opts->{g} // 1; #/

    my $cmd = $CONF->{TOPHAT};
    die "ERROR: TOPHAT is not in the config." unless ($cmd);

    my $args = [
        ["-o", ".", 0],
        ["-g", $g, 0],
        ["-p", '32', 0],
        ["", $name, 1],
        ["", $fastq, 1]
    ];

    my $inputs = [
        $fasta,
        $fastq,
    ];

    # add gff file if genome has annotations
    unshift @$args, ["-G", $gff, 1] if $gff;
    unshift @$inputs, $gff if $gff;
    
    # add index files
    push @$inputs, @index_files;

    return (
        cmd => $cmd,
        script => undef,
        options => {
            "allow-zero-length" => JSON::false,
        },
        args => $args,
        inputs => $inputs,
        outputs => [
            catfile($staging_dir, "accepted_hits.bam")
        ],
        description => "Running tophat..."
    );
}

sub create_gsnap_workflow {
    my %opts = @_;

    # Required arguments
    my $gid = $opts{gid};
    my $fasta = $opts{fasta};
    my $fastq = $opts{fastq};
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $alignment_params = $opts{alignment_params} // {}; #/

    # Generate index
    my %gmap = create_gmap_index_job($gid, $fasta);

    # Generate sam file
    my %gsnap = create_gsnap_job({
        fastq => $fastq,
        gmap => @{@{$gmap{outputs}}[0]}[0],
        staging_dir => $staging_dir,
        alignment_params => $alignment_params,
    });

    # Generate and sort bam
    my %bam = create_samtools_bam_job(@{$gsnap{outputs}}[0], $staging_dir);
    my %sorted_bam = create_samtools_sort_job(@{$bam{outputs}}[0], $staging_dir);

    # Return the bam output name and jobs required
    return @{$sorted_bam{outputs}}[0], (
        \%gmap,
        \%gsnap,
        \%bam,
        \%sorted_bam
    );
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
    my $name = to_filename($samfile);
    my $cmd = $CONF->{SAMTOOLS};

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["view", '', 0],
            ["-bS", $samfile, 1],
            [">", $name . ".bam", 0]
        ],
        inputs => [
            $samfile
        ],
        outputs => [
            catfile($staging_dir, $name . ".bam")
        ],
        description => "Generating bam file..."
    );
}

sub create_samtools_sort_job {
    my $bam = shift;
    my $staging_dir = shift;
    my $name = to_filename($bam);
    my $cmd = $CONF->{SAMTOOLS};

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["sort", '', 0],
            ["", $bam, 1],
            ["", $name . "-sorted", 1]
        ],
        inputs => [
            $bam
        ],
        outputs => [
            catfile($staging_dir, $name . "-sorted.bam")
        ],
        description => "Sorting bam file..."
    );
}

sub create_gsnap_job {
    my $opts = shift;

    # Required arguments
    my $fastq = $opts->{fastq};
    my $gmap = $opts->{gmap};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $aligner = $opts->{aligner};
    my $gapmode = $aligner->{gap} // "none"; #/
    my $Q = $aligner->{Q} // 1; #/
    my $n = $aligner->{n} // 5; #/
    my $nofail = $aligner->{nofail} // 1; #/

    my $name = basename($gmap);
    my $cmd = $CONF->{GSNAP};

    my $args = [
        ["-D", ".", 0],
        ["-d", $name, 0],
        ["--nthreads=32", '', 0],
        ["-n", $n, 0],
        ["--format=sam", '', 0],
        ["--gmap-mode=$gapmode", "", 1],
        ["--batch=5", $fastq, 0],
    ];

    push $args, ["-Q", "", 0] if $Q;
    push $args, ["--nofails", "", 1] if $nofail;

    return (
        cmd => qq[$cmd > $name.sam],
        script => undef,
        options => {
            "allow-zero-length" => JSON::false,
        },
        args => $args,
        inputs => [
            $fastq,
            [$gmap, 1]
        ],
        outputs => [
            catfile($staging_dir, $name . ".sam")
        ],
        description => "Running gsnap..."
    );
}

sub create_alignment_workflow { #TODO make this into a pipeline module
    my %opts = @_;

    # Required arguments
    my $user = $opts{user};
    my $wid = $opts{wid}; # workflow ID
    my $fastq = $opts{input_file}; # input fastq file
    my $genome = $opts{genome}; # genome object
    my $staging_dir = $opts{staging_dir};
    my $result_dir = $opts{result_dir};
    my $metadata = $opts{metadata};
    my $options = $opts{options}; # Options for cutadapt and aligner
    my $alignment_type = $options->{alignment_params}->{tool}; # "tophat" or "gsnap"
    #print STDERR "CommonTasks::create_alignment_workflow ", Dumper $options, "\n";

    my $gid = $genome->id;
    my $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;
    my @tasks;

    # Validate the fastq file
    push @tasks, create_validate_fastq_job($fastq);

    # Reheader the fasta file
    my $fasta = get_genome_file($genome->id);
    my $reheader_fasta = to_filename($fasta) . ".reheader.faa";
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta,
        reheader_fasta => $reheader_fasta,
        cache_dir => $FASTA_CACHE_DIR
    );

    # Index the fasta file
    push @tasks, create_fasta_index_job(
        fasta => catfile($FASTA_CACHE_DIR, $reheader_fasta),
        cache_dir => $FASTA_CACHE_DIR
    );

    # Trim the fastq file
    my $trimmed = create_cutadapt_job(
        fastq => $fastq,
        validated => "$fastq.validated",
        staging_dir => $staging_dir,
        cutadapt => $options->{alignment_params}->{cutadapt_params},
    );
    my $trimmed_fastq = @{$trimmed->{outputs}}[0];
    push @tasks, $trimmed;

    # Generate gff if genome annotated
    my $gff_file;
    if ( $genome->has_gene_features ) {
        my $gff = create_gff_generation_job(
            gid => $gid,
            organism_name => $genome->organism->name,
            #validated => "$fastq.validated"
        );
        $gff_file = @{$gff->{outputs}}[0];
        push @tasks, $gff;
    }

    # Add aligner workflow
    my ($bam, @alignment_tasks);
    if ($alignment_type eq 'tophat') {
        ($bam, @alignment_tasks) = create_tophat_workflow(
            gid => $gid,
            fasta => catfile($FASTA_CACHE_DIR, $reheader_fasta),
            fastq => $trimmed_fastq,
            gff => $gff_file,
            staging_dir => $staging_dir,
            aligner_params => $options->{alignment_params},
        );
    }
    elsif ($alignment_type eq 'gsnap') {
        ($bam, @alignment_tasks) = create_gsnap_workflow(
            gid => $gid,
            fasta => $reheader_fasta,
            fastq => $trimmed_fastq,
            staging_dir => $staging_dir,
            alignment_params => $options->{alignment_params},
        );
    }
    else {
        print STDERR "Error: unrecognized alignment '$alignment'\n";
        return;
    }
    push @tasks, @alignment_tasks;

    # Load alignment
    push @tasks, create_load_bam_job(
        user => $user,
        metadata => $metadata,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        annotations => '', #FIXME 12/12/14
        wid => $wid,
        gid => $gid,
        bam_file => $bam
    );

    return (\@tasks, $bam, $reheader_fasta, $gff_file);
}

1;
