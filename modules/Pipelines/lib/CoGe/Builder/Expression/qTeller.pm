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

use CoGe::Accessory::Jex;
use CoGe::Accessory::TDS qw(read);
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Workflow;
use CoGe::Core::Storage qw(get_genome_file get_experiment_files get_workflow_paths);
use CoGe::Builder::CommonTasks;

our $CONF = CoGe::Accessory::Web::get_defaults();

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
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
    my %opts = @_;
    my $genome = $opts{genome};
    my $user = $opts{user};
    my $input_file = $opts{input_file}; # path to bam file
    my $metadata = $opts{metadata};
    my $staging_dir = $opts{staging_dir};
    my $result_dir = $opts{result_dir};
    my $wid = $opts{wid};
    my $options = $opts{options};
    #print STDERR "qTeller::build ", Dumper $metadata, "\n";

    my $gid = $genome->id;
    my $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;

    # Check if genome has annotations
    my $isAnnotated = $genome->has_gene_features;

    # Set metadata for the pipeline being used
    my $annotations = ""; #FIXME generate_metadata($alignment, $isAnnotated);

    my @jobs;

    # Filter the fasta file (clean up headers)
    my $fasta = get_genome_file($gid);
    my $reheader_fasta = to_filename($fasta) . ".reheader.fasta";
    push @jobs, create_fasta_reheader_job( fasta => $fasta, reheader_fasta => $reheader_fasta, cache_dir => $FASTA_CACHE_DIR );
    
    # Generate gff if genome annotated
    my $gff_file;
    if ($isAnnotated) {
        my $gff = create_gff_generation_job(gid => $gid, organism_name => $genome->organism->name);
        $gff_file = @{$gff->{outputs}}[0];
        push @jobs, $gff;
    }

    # Generate bed file
    my $bed = create_bed_file_job(
        bam => $input_file,
        staging_dir => $staging_dir,
        seq => $options->{expression_params},
    );
    push @jobs, $bed;

    # Filter bed file
    my $filtered_bed = create_filter_bed_file_job(@{$bed->{outputs}}[0], $staging_dir);
    push @jobs, $filtered_bed;

    # Check for annotations required by cufflinks
    my $include_csv;
    if ($isAnnotated) {
        $include_csv = 1;

        # Run cufflinks
        my $cuff = create_cufflinks_job($gff_file, catfile($FASTA_CACHE_DIR, $reheader_fasta), $input_file, $staging_dir);
        push @jobs, $cuff;

        # Convert final output into csv
        my $parse_cuff = create_parse_cufflinks_job(@{$cuff->{outputs}}[0], $staging_dir);
        push @jobs, $parse_cuff;

        # Load csv experiment
        my $load_csv = create_load_csv_job($metadata, $gid, @{$parse_cuff->{outputs}}[0], $user, $annotations, $staging_dir, $wid);
        push @jobs, $load_csv;
    }

    # Load bam experiment
    my $load_bam = create_load_bam_job($metadata, $gid, $input_file, $user, $annotations, $staging_dir, $wid);
    push @jobs, $load_bam;

    # Load bed experiment
    my $load_bed = create_load_bed_job($metadata, $gid, @{$filtered_bed->{outputs}}[0], $user, $annotations, $staging_dir, $wid);
    push @jobs, $load_bed;

    # Create notebook
    my $notebook = create_notebook_job($metadata, $include_csv, 1, 1, $user, $annotations, $staging_dir, $result_dir);
    push @jobs, $notebook;

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
    }
    else {
        push @annotations, (
            qq{note|gmap_build},
            qq{note|gsnap -n 5 --format=sam -Q --gmap-mode=none --nofails},
            qq{note|samtools mpileup -D -Q 20}
        );
    }

    push @annotations, qq{note|cufflinks} if $annotated;

    return join(';', @annotations);
}

# mdb removed 12/12/14 - replaced by CoGeX::Genome::has_gene_features
#sub has_annotations {
#    my ($gid, $db) = @_;
#
#    #FIXME Remove hardcoded values
#    my $count = $db->resultset('Feature')->count(
#        {
#            feature_type_id => [ 1, 2, 3 ],
#            genome_id       => $gid
#        },
#        { join => [ { dataset => 'dataset_connectors' } ], }
#    );
#
#    return $count > 0;
#}

sub create_cufflinks_job {
    print STDERR "cufflinks ", Dumper \@_, "\n";
    my ($gff, $fasta, $bam, $staging_dir) = @_;
    my $cmd = $CONF->{CUFFLINKS};
    die "ERROR: CUFFLINKS is not in the config." unless ($cmd);

    return {
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
    };
}

sub create_bed_file_job {
    my %opts = @_;

    # Required arguments
    my $bam = $opts{bam};
    my $staging_dir = $opts{staging_dir};

    # Optional arguments
    my $seq = $opts{seq} // {};
    my $Q = $seq->{Q} // 20;

    my $name = to_filename($bam);
    my $cmd = $CONF->{SAMTOOLS};
    my $PILE_TO_BED = catfile($CONF->{SCRIPTDIR}, "pileup_to_bed.pl");
    die "ERROR: SAMTOOLS is not in the config." unless ($cmd);
    die "ERROR: SCRIPTDIR not specified in config" unless $PILE_TO_BED;

    return {
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
    };
}

sub create_filter_bed_file_job {
    my $bed = shift;
    my $staging_dir = shift;
    my $name = to_filename($bed);
    my $cmd = $CONF->{SAMTOOLS};
    my $NORMALIZE_BED = catfile($CONF->{SCRIPTDIR}, "normalize_bed.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $NORMALIZE_BED;

    return {
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
    };
}

sub create_parse_cufflinks_job {
    my $cufflinks = shift;
    my $staging_dir = shift;

    my $name = to_filename($cufflinks);

    my $cmd = $CONF->{PYTHON};
    my $script = $CONF->{PARSE_CUFFLINKS};
    die "ERROR: PYTHON not in the config." unless ($cmd);
    die "ERROR: PARSE_CUFFLINKS is not in the config." unless ($script);

    return {
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
    };
}

sub create_load_csv_job {
    my ($md, $gid, $csv, $user, $annotations, $staging_dir, $wid) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return {
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
            ['-source_name', '"'.$md->{source}.'"', 0],
            ['-types', qq{"Expression"}, 0],
            ['-annotations', qq{"$annotations"}, 0],
            ['-staging_dir', "./csv", 0],
            ['-file_type', "csv", 0],
            ['-data_file', $csv, 0],
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
    };
}

sub create_load_bam_job {
    my ($md, $gid, $bam, $user, $annotations, $staging_dir, $wid) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return {
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
            ['-source_name', '"'.$md->{source}.'"', 0],
            ['-types', qq{"RNAseq;BAM"}, 0],
            ['-annotations', qq{"$annotations"}, 0],
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
    };
}

sub create_load_bed_job {
    my ($md, $gid, $bed, $user, $annotations, $staging_dir, $wid) = @_;
    my $cmd = catfile($CONF->{SCRIPTDIR}, "load_experiment.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return {
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
            ['-source_name', '"'.$md->{source}.'"', 0],
            ['-types', qq{"Expression"}, 0],
            ['-annotations', qq{"$annotations"}, 0],
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
    };
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
        ['-annotations', qq{"$annotations"}, 0],
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

    return {
        cmd => $cmd,
        script => undef,
        args => $args,
        inputs => $inputs,
        outputs => [
            catfile($result_dir, '1')
        ],
        description => "Creating notebook..."
    };
}

1;
