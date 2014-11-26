package CoGe::Builder::Expression::qTeller;

use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Basename qw(fileparse basename dirname);
use File::Path qw(mkpath);
use File::Spec::Functions qw(catdir catfile);
use JSON qw(decode_json);
use URI::Escape::JavaScript qw(unescape);

use CoGe::Accessory::TDS qw(read);
use CoGe::Accessory::Workflow;
use CoGe::Accessory::Jex;
use CoGe::Core::Storage qw(get_genome_file get_experiment_files get_workflow_paths);
use CoGe::Accessory::Web qw(get_defaults);

our $CONF = CoGe::Accessory::Web::get_defaults();

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT $CACHE @EXPORT_OK);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw( run );
    @EXPORT_OK = qw(build);
}

sub run {
    my $opts = shift;
    my $user = $opts->{user};

    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Running the qTeller pipeline', init => 1 );
    my $wid = $workflow->id;

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    my $options = {
        result_dir => $result_dir,
        staging_dir => $staging_dir,
        wid  => $wid,
        %{$opts},
    };

    my @jobs = build($options);

    for my $job (@jobs) {
        $workflow->add_job(%{$job});
    }

    #say STDERR "WORKFLOW DUMP\n" . Dumper($workflow) if $DEBUG;
    #say STDERR "JOB NOT SCHEDULED TEST MODE" and exit(0) if $test;

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub build {
    my $opts = shift;
    my $db = $opts->{db};
    my $genome = $opts->{genome};
    my $user = $opts->{user};
    my $files = $opts->{files};
    my $metadata = $opts->{metadata};
    my $alignment = $opts->{alignment_type};
    my $staging_dir = $opts->{staging_dir};
    my $result_dir = $opts->{result_dir};
    my $wid = $opts->{wid};
    my $options = $opts->{options};

    my $gid = $genome->id;

    $CACHE = $CONF->{CACHEDIR};
    die "ERROR: CACHEDIR not specified in config" unless $CACHE;
    mkpath($CACHE, 0, 0777) unless -r $CACHE;

    # Set default alignment program if alignment is not set or incorrect
    $alignment = "gsnap" unless ($alignment and $alignment =~ /tophat|gsnap/);

    # Check if genome has annotations
    my $annotated = has_annotations($genome->id, $db);

    # Set metadata for the pipeline being used
    my $annotations = generate_metadata($alignment, $annotated);

    my $fasta = get_genome_file($gid);
    my $input_file = shift(@$files);
    my $file_type = $input_file->{type};
    my $fastq = $input_file->{path};

    my (@jobs, @steps, $bam, $include_csv);

    # Filter the fasta file (clean up headers)
    my %filter = create_fasta_filter_job($gid, $fasta, "$fastq.validated");
    my $filtered_fasta = @{$filter{outputs}}[0];
    push @jobs, \%filter;

    # Validate the fastq file
    # XXX: Add job dependencies
    my %validate = create_validate_fastq_job($fastq);
    push @jobs, \%validate;

    # Cleanup fastq
    my %trimmed = create_cutadapt_job({
        fastq => $fastq,
        validated => "$fastq.validated",
        staging_dir => $staging_dir,
        cutadapt => $options->{cutadapt},
    });

    my $trimmed_fastq = @{$trimmed{outputs}}[0];
    push @jobs, \%trimmed;

    my $gff_file;
    # Generate gff if genome annotated
    if ($annotated) {
        my %gff = create_gff_generation_job($genome, "$fastq.validated");
        $gff_file = @{$gff{outputs}}[0];
        push @jobs, \%gff;
    }

    if ($alignment eq "tophat") {
        ($bam, @steps) = tophat_pipeline({
            gid => $gid,
            fasta => $filtered_fasta,
            fastq => $trimmed_fastq,
            gff => $gff_file,
            staging_dir => $staging_dir,
            aligner => $options->{aligner}
        });
    }
    else {
        ($bam, @steps) = gsnap_pipeline({
            gid => $gid,
            fasta => $filtered_fasta,
            fastq => $trimmed_fastq,
            staging_dir => $staging_dir,
            aligner => $options->{aligner}
        });
    }

    # Join alignment pipeline
    @jobs = (@jobs, @steps);

    # Generate bed file
    my %bed = create_bed_file_job({
        bam => $bam,
        staging_dir => $staging_dir,
        seq => $options->{seq},
    });

    push @jobs, \%bed;

    # Filter bed file
    my %filtered_bed = create_filter_bed_file_job(@{$bed{outputs}}[0], $staging_dir);
    push @jobs, \%filtered_bed;

    # Check for annotations required by cufflinks
    if ($annotated) {
        $include_csv = 1;

        # Run cufflinks
        my %cuff = create_cufflinks_job($gff_file, $filtered_fasta, $bam, $staging_dir);
        push @jobs, \%cuff;

        # Convert final output into csv
        my %parse_cuff = create_parse_cufflinks_job(@{$cuff{outputs}}[0], $staging_dir);
        push @jobs, \%parse_cuff;

        # Load csv experiment
        my %load_csv = create_load_csv_job($metadata, $gid, @{$parse_cuff{outputs}}[0], $user, $annotations, $staging_dir, $wid);
        push @jobs, \%load_csv;
    }

    # Load bam experiment
    my %load_bam = create_load_bam_job($metadata, $gid, $bam, $user, $annotations, $staging_dir, $wid);
    push @jobs, \%load_bam;

    # Load bed experiment
    my %load_bed = create_load_bed_job($metadata, $gid, @{$filtered_bed{outputs}}[0], $user, $annotations, $staging_dir, $wid);
    push @jobs, \%load_bed;

    # Create notebook
    my %notebook = create_notebook_job($metadata, $include_csv, 1, 1, $user, $annotations, $staging_dir, $result_dir);
    push @jobs, \%notebook;

    return wantarray ? @jobs : \@jobs;
}

sub generate_metadata {
    my ($pipeline, $annotated) = @_;

    my @annotations = (
        qq{http://genomevolution.org/wiki/index.php/Expression_Analysis_Pipeline||note|Generated by CoGe's Expression Analysis Pipeline},
        qq{note|cutadapt -q 25 -m 17},
    );

    if ($pipeline eq "tophat") {
        push @annotations, (
            qq{note|bowtie2_build},
            qq{note|tophat -g 1}
        );
    } else {
        push @annotations, (
            qq{note|gmap_build},
            qq{note|gsnap -n 5 --format=sam -Q --gmap-mode=none --nofails},
            qq{note|samtools mpileup -D -Q 20}
        );
    }

    push @annotations, qq{note|cufflinks} if $annotated;

    return '"' . join(';', @annotations) . '"';
}

sub to_filename {
    my ($name, undef, undef) = fileparse(shift, qr/\.[^.]*/);
    return $name;
}

sub sanitize_organism_name {
    my $org = shift;

    $org =~ s/\///g;
    $org =~ s/\s+/_/g;
    $org =~ s/\(//g;
    $org =~ s/\)//g;
    $org =~ s/://g;
    $org =~ s/;//g;
    $org =~ s/#/_/g;
    $org =~ s/'//g;
    $org =~ s/"//g;

    return $org;
}

sub has_annotations {
    my ($gid, $db) = @_;

    #FIXME Remove hardcoded values
    my $count = $db->resultset('Feature')->count(
        {
            feature_type_id => [ 1, 2, 3 ],
            genome_id       => $gid
        },
        { join => [ { dataset => 'dataset_connectors' } ], }
    );

    return $count > 0;
}

sub create_validate_fastq_job {
    my $fastq = shift;

    my $cmd = catfile($CONF->{SCRIPTDIR}, "validate_fastq.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
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
    );
}

sub create_fasta_filter_job {
    my $gid = shift;
    my $fasta = shift;
    my $validated = shift;
    my $name = to_filename($fasta);
    my $cmd = catfile($CONF->{SCRIPTDIR}, "fasta_reheader.pl");
    my $FASTA_CACHE_DIR = catdir($CACHE, $gid, "fasta");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fasta, 1],
            ["", $name . ".filtered.fa", 0]
        ],
        inputs => [
            $fasta,
            $validated
        ],
        outputs => [
            catfile($FASTA_CACHE_DIR, $name . ".filtered.fa")
        ],
        description => "Filtering genome sequence..."
    );
}

sub create_gff_generation_job {
    my $genome = shift;
    my $validated = shift;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "coge_gff.pl");
    my $org_name = sanitize_organism_name($genome->organism->name);
    my $name = "$org_name-1-name-0-0-id-" . $genome->id . "-1.gff";

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-f', $name, 0],
            ['-staging_dir', '.', 0],
            ['-gid', $genome->id, 0],
            ['-upa', 1, 0],
            ['-id_type', "name", 0],
            ['-cds', 0, 0],
            ['-annos', 0, 0],
            ['-nu', 1, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1],
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $validated
        ],
        outputs => [
            catdir($CACHE, $genome->id, "gff", $name)
        ],
        description => "Generating genome annotations gff..."
    );
}

sub create_cutadapt_job {
    my $opts = shift;

    # Required arguments
    my $fastq = $opts->{fastq};
    my $validated = $opts->{validated};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $cutadapt= $opts->{cutadapt} // {};
    my $q = $cutadapt->{q} // 25;
    my $quality = $cutadapt->{quality} // 64;
    my $m = $cutadapt->{m} // 17;

    my $name = to_filename($fastq);
    my $cmd = $CONF->{CUTADAPT};
    die "ERROR: CUTADAPT is not in the config." unless ($cmd);

    return (
        cmd => qq[$cmd > /dev/null],
        script => undef,
        args => [
            ['-q', $q, 0],
            ["--quality-base=$quality", '', 0],
            ['-m', $m, 0],
            ['', $fastq, 1],
            ['-o', $name . '.trimmed.fastq', 1],
        ],
        inputs => [
            $fastq,
            $validated
        ],
        outputs => [
            catfile($staging_dir, $name . '.trimmed.fastq')
        ],
        description => "Running cutadapt..."
    );
}

sub create_cufflinks_job {
    my ($gff, $fasta, $bam, $staging_dir) = @_;
    my $cmd = $CONF->{CUFFLINKS};
    die "ERROR: CUFFLINKS is not in the config." unless ($cmd);

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-q', '', 0], # suppress output other than warning/error messages
            ['-u', '', 0],
            ['-b', $fasta, 1],
            ['-p', 24, 0],
            ['', $bam, 1]
        ],
        inputs => [
            $gff,
            $bam,
            $fasta
        ],
        outputs => [
            catfile($staging_dir, "genes.fpkm_tracking")
        ],
        description => "Running cufflinks..."
    );
}

sub create_bed_file_job {
    my $opts = shift;

    # Required arguments
    my $bam = $opts->{bam};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $seq = $opts->{seq} // {};
    my $Q = $seq->{Q} // 20;

    my $name = to_filename($bam);
    my $cmd = $CONF->{SAMTOOLS};
    my $PILE_TO_BED = catfile($CONF->{SCRIPTDIR}, "pileup_to_bed.pl");
    die "ERROR: SAMTOOLS is not in the config." unless ($cmd);
    die "ERROR: SCRIPTDIR not specified in config" unless $PILE_TO_BED;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['mpileup', '', 0],
            ['-D', '', 0],
            ['-Q', $Q, 0],
            ['', $bam, 1],
            ['|', 'perl', 0],
            [$PILE_TO_BED, '', 0],
            ['>', $name . ".bed",  0]
        ],
        inputs => [
            $bam,
        ],
        outputs => [
            catfile($staging_dir, $name . ".bed")
        ],
        description => "Generating read depth..."
    );
}

sub create_filter_bed_file_job {
    my $bed = shift;
    my $staging_dir = shift;
    my $name = to_filename($bed);
    my $cmd = $CONF->{SAMTOOLS};
    my $NORMALIZE_BED = catfile($CONF->{SCRIPTDIR}, "normalize_bed.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $NORMALIZE_BED;

    return (
        cmd => "perl",
        script => undef,
        args => [
            [$NORMALIZE_BED, $bed, 0],
            ['>', $name . '.normalized.bed', 0]
        ],
        inputs => [
            $bed,
        ],
        outputs => [
            catfile($staging_dir, $name . ".normalized.bed")
        ],
        description => "Normalizing read depth..."
    );
}

sub create_parse_cufflinks_job {
    my $cufflinks = shift;
    my $staging_dir = shift;
    my $name = to_filename($cufflinks);
    my $cmd = $CONF->{PYTHON};
    my $script = $CONF->{PARSE_CUFFLINKS};
    die "ERROR: PYTHON not in the config." unless ($cmd);
    die "ERROR: PARSE_CUFFLINKS is not in the config." unless ($script);

    return (
        cmd => "$cmd $script",
        script => undef,
        args => [
            ["", $cufflinks, 0],
            ["", $name . ".csv", 0]
        ],
        inputs => [
            $cufflinks
        ],
        outputs => [
            catfile($staging_dir, $name . ".csv")
        ],
        description => "Processing cufflinks output ..."
    );
}

sub create_load_csv_job {
    my ($md, $gid, $csv, $user, $annotations, $staging_dir, $wid) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', '"'.$md->{name}.' (FPKM)'.'"', 0],
            ['-desc', qq{"Transcript expression measurements"}, 0],
            ['-version', '"'.$md->{version}.'"', 0],
            ['-restricted', $md->{restricted}, 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', '"'.$md->{source_name}.'"', 0],
            ['-types', qq{"Expression"}, 0],
            ['-annotations', $annotations, 0],
            ['-staging_dir', "./csv", 0],
            ['-file_type', "csv", 0],
            ['-data_file', "$csv", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $csv
        ],
        outputs => [
            [catdir($staging_dir, "csv"), 1],
            catfile($staging_dir, "csv/log.done"),
        ],
        description => "Loading expression data..."
    );
}

sub create_load_bam_job {
    my ($md, $gid, $bam, $user, $annotations, $staging_dir, $wid) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', '"'.$md->{name}." (alignment)".'"', 0],
            ['-desc', qq{"Mapped reads"}, 0],
            ['-version', '"'.$md->{version}.'"', 0],
            ['-restricted', $md->{restricted}, 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', '"'.$md->{source_name}.'"', 0],
            ['-types', qq{"RNAseq;BAM"}, 0],
            ['-annotations', $annotations, 0],
            ['-staging_dir', "./bam", 0],
            ['-file_type', "bam", 0],
            ['-data_file', "$bam", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $bam,
        ],
        outputs => [
            [catdir($staging_dir, "bam"), 1],
            catfile($staging_dir, "bam/log.done"),
        ],
        description => "Loading mapped reads..."
    );
}

sub create_load_bed_job {
    my ($md, $gid, $bed, $user, $annotations, $staging_dir, $wid) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', '"'.$md->{name}." (read depth)".'"', 0],
            ['-desc', qq{"Read depth per position"}, 0],
            ['-version', '"'.$md->{version}.'"', 0],
            ['-restricted', $md->{restricted}, 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', '"'.$md->{source_name}.'"', 0],
            ['-types', qq{"Expression"}, 0],
            ['-annotations', $annotations, 0],
            ['-staging_dir', "./bed", 0],
            ['-file_type', "bed", 0],
            ['-data_file', "$bed", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $bed,
        ],
        outputs => [
            [catdir($staging_dir, "bed"), 1],
            catfile($staging_dir, "bed/log.done"),
        ],
        description => "Loading read depth..."
    );
}

sub create_notebook_job {
    my ($md, $csv, $bam, $bed, $user, $annotations, $staging_dir, $result_dir) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "create_notebook.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    my $args = [
        ['-uid', $user->id, 0],
        ['-name', '"'.$md->{name}.'"', 0],
        ['-desc', '"'.$md->{description}.'"', 0],
        ['-page', '""', 0],
        ['-type', 2, 0],
        ['-restricted', $md->{restricted}, 0],
        ['-annotations', $annotations, 0],
        ['-config', $CONF->{_CONFIG_PATH}, 1],
        ['-log', catfile($staging_dir, "log.txt"), 0],
        ['-result_dir', $result_dir, 0]
    ];

    my $inputs = [$CONF->{_CONFIG_PATH}];

    if ($bed) {
        push @$args, ['', "bed/log.txt", 0];
        push @$inputs, [catdir($staging_dir, "bed"), 1];
        push @$inputs, catfile($staging_dir, "bed/log.done"),
    }

    if ($bam) {
        push @$args, ['', "bam/log.txt", 0];
        push @$inputs, [catdir($staging_dir, "bam"), 1];
        push @$inputs, catfile($staging_dir, "bam/log.done");
    }

    if ($csv) {
        push @$args, ['', "csv/log.txt", 0];
        push @$inputs, [catdir($staging_dir, "csv"), 1];
        push @$inputs, catfile($staging_dir, "csv/log.done");
    }

    return (
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => $inputs,
        outputs => [
            catfile($result_dir, '1')
        ],
        description => "Creating notebook..."
    );
}

#
# GSNAP PIPELINE AND JOBS
#
sub gsnap_pipeline {
    my $opts = shift;

    # Required arguments
    my $gid = $opts->{gid};
    my $fasta = $opts->{fasta};
    my $fastq = $opts->{fastq};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $aligner = $opts->{aligner} // {};

    # Generate index
    my %gmap = create_gmap_index_job($gid, $fasta);

    # Generate sam file
    my %gsnap = create_gsnap_job({
        fastq => $fastq,
        gmap => @{@{$gmap{outputs}}[0]}[0],
        staging_dir => $staging_dir,
        aligner => $aligner,
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
    my $GMAP_CACHE_DIR = catdir($CACHE, $gid, "gmap_index");
    die "ERROR: GMAP_BUILD is not in the config." unless ($cmd);

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

sub create_gsnap_job {
    my $opts = shift;

    # Required arguments
    my $fastq = $opts->{fastq};
    my $gmap = $opts->{gmap};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $aligner = $opts->{aligner};
    my $gapmode = $aligner->{gap} // "none";
    my $Q = $aligner->{Q} // 1;
    my $n = $aligner->{n} // 5;
    my $nofail = $aligner->{nofail} // 1;

    my $name = basename($gmap);
    my $cmd = $CONF->{GSNAP};
    die "ERROR: GSNAP is not in the config." unless ($cmd);

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

sub create_samtools_bam_job {
    my $samfile = shift;
    my $staging_dir = shift;
    my $name = to_filename($samfile);
    my $cmd = $CONF->{SAMTOOLS};
    die "ERROR: SAMTOOLS is not in the config." unless ($cmd);

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
    die "ERROR: SAMTOOLS is not in the config." unless ($cmd);

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

#
# TopHat pipeline
#
sub tophat_pipeline {
    my $opts = shift;

    # Required arguments
    my $gid = $opts->{gid};
    my $fasta = $opts->{fasta};
    my $fastq = $opts->{fastq};
    my $gff = $opts->{gff};
    my $staging_dir = $opts->{staging_dir};

    # Optional arguments
    my $aligner = $opts->{aligner} // {};

    my $BOWTIE_CACHE_DIR = catdir($CACHE, $gid, "bowtie_index");

    my ($index, %bowtie) = create_bowtie_index_job($gid, $fasta);
    my %tophat = create_tophat_job({
        staging_dir => $staging_dir,
        fastq => $fastq,
        fasta => $fasta,
        gff   => $gff,
        index_name => $index,
        index_files => ($bowtie{outputs}),
        g =>  $aligner->{g},
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
    my $BOWTIE_CACHE_DIR = catdir($CACHE, $gid, "bowtie_index");
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
    my $g = $opts->{g} // 1;

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

1;
