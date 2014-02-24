#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Basename qw(fileparse basename);
use File::Path qw(mkpath);
use File::Spec qw(catdir);
use Getopt::Long qw(GetOptions);
use JSON qw(decode_json);
use URI::Escape::JavaScript qw(unescape);

use CoGe::Accessory::Workflow;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Storage qw(get_genome_file get_experiment_files);
use CoGe::Accessory::Web qw(get_defaults get_job schedule_job);

our ($DESC, $YERBA, $LOG, $DEBUG, $P, $db, $host, $port, $user, $pass,
     $name, $description, $version, $restricted, $source_name, $files,
     $CACHE, $test, $config, $gid, $jobid, $userid, $staging_dir);

$DESC = "Running the qTeller pipeline";

GetOptions(
    "debug|d=s"         => \$DEBUG,       # Dumps the workflow hash
    "test|t=s"          => \$test,        # Skip workflow submission
    "data_file|df=s"    => \$files,       # Input files to be processed
    "gid=s"             => \$gid,         # Reference genome
    "jobid|jid=s"       => \$jobid,       # Reference job
    "userid|uid=s"      => \$userid,      # User loading the experiment

    # General configuratino options
    "staging_dir|dir=s" => \$staging_dir,
    "config|cfg=s"      => \$config,

    # Load experiment options
    "name|n=s"          => \$name,
    "description|d=s"   => \$description,
    "version|v=s"       => \$version,
    "restricted|r=s"    => \$restricted,
    "source_name|s=s"   => \$source_name
);

$| = 1;

# Enforce mandatory options
die "ERROR: genome id not specified use gid" unless $gid;
die "ERROR: job id not specified use jobid" unless $jobid;
die "ERROR: user id not specified use userid" unless $userid;
die "ERROR: staging directory not specified use staging_dir or s" unless $staging_dir;
die "ERROR: config not specified use config or cfg" unless $config;
die "ERROR: no experiment files found" unless $files;

sub setup {
    $P    = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};
    $YERBA = CoGe::Accessory::Jex->new( host => $P->{JOBSERVER}, port => $P->{JOBPORT} );
    $LOG = File::Spec->catdir(($staging_dir, "qteller.txt"));

    $CACHE = $P->{CACHEDIR};
    die "ERROR: CACHEDIR not specified in config" unless $CACHE;

    mkpath($CACHE, 0, 0777) unless -r $CACHE;
    mkpath($staging_dir, 0, 0777) unless -r $staging_dir;
    my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";

    return CoGeX->connect( $connstr, $user, $pass );
}

sub main {
    my $coge = setup();
    die "ERROR: couldn't connect to the database" unless $coge;

    my $user = $coge->resultset("User")->find($userid);
    die "ERROR: user could not be found for id: $userid" unless $user;

    my $job = $coge->resultset("Job")->find($jobid);
    die "ERROR: the job could not be fetched for $jobid" unless $job;

    my $genome = $coge->resultset("Genome")->find($gid);
    die "ERROR: a genome with $gid not found" unless $genome;

    my $fasta = get_genome_file($gid);

    #XXX: files should allow for comma seperated list of files
    my $fastq = $files;
    die "ERROR: no experiment files found" unless $fastq;

    my $workflow = $YERBA->create_workflow(
        id => $job->id,
        name => $DESC,
        logfile => $LOG
    );


    # Validate the fastq file
    # XXX: Add job dependencies
    my %validate = create_validate_fastq_job($fastq);

    # Filter the fasta file (clean up headers)
    my %filtered_fasta = create_fasta_filter_job($fasta, "$fastq.validated");

    # Generate gff
    my %gff = create_gff_generation_job($genome, "$fastq.validated");

    # Cleanup fastq
    my %trimmed = create_cutadapt_job($fastq, "$fastq.validated");

    # Generate index
    my %gmap = create_gmap_index_job(@{$filtered_fasta{outputs}}[0]);

    # Generate sam file
    my %gsnap = create_gsnap_job(@{$trimmed{outputs}}[0], @{@{$gmap{outputs}}[0]}[0]);

    # Generate and sort bam
    my %bam = create_samtools_bam_job(@{$gsnap{outputs}}[0]);
    my %sorted_bam = create_samtools_sort_job($bam{outputs}[0]);

    # Run cufflinks
    my %cuff = create_cufflinks_job(@{$gff{outputs}}[0],
        @{$filtered_fasta{outputs}}[0],
        @{$sorted_bam{outputs}}[0]);

    # Convert final output into csv
    my %parse_cuff = create_parse_cufflinks_job(@{$cuff{outputs}}[0]);

    # Load csv experiment
    my %load_csv = create_load_csv_job(@{$parse_cuff{outputs}}[0], $user);

    # Load bam experiment
    my %load_bam = create_load_bam_job(@{$parse_cuff{outputs}}[0],
        @{$sorted_bam{outputs}}[0],
        $user);

    my %create_notebook = create_notebook_job(@{$load_csv{outputs}}[0],
        @{$load_bam{outputs}}[0],
        $user);

    $workflow->add_job(%validate);
    $workflow->add_job(%filtered_fasta);
    $workflow->add_job(%gff);
    $workflow->add_job(%trimmed);
    $workflow->add_job(%gmap);
    $workflow->add_job(%gsnap);
    $workflow->add_job(%bam);
    $workflow->add_job(%sorted_bam);
    $workflow->add_job(%cuff);
    $workflow->add_job(%parse_cuff);
    $workflow->add_job(%load_csv);
    $workflow->add_job(%load_bam);
    $workflow->add_job(%create_notebook);

    say STDERR "WORKFLOW DUMP\n" . Dumper($workflow) if $DEBUG;
    say STDERR "JOB NOT SCHEDULED TEST MODE" and exit(0) if $test;

    my $result = $YERBA->submit_workflow($workflow);
    # check if the schedule was successful
    my $status = decode_json($result);
    exit(1) if defined($status->{error}) and lc($status->{error}) eq "error";
    CoGe::Accessory::Web::schedule_job(job => $job);
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

sub create_validate_fastq_job {
    my $fastq = shift;

    my $cmd = File::Spec->catdir(($P->{SCRIPTDIR}, "validate_fastq.pl"));
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
    my $fasta = shift;
    my $validated = shift;
    my $name = to_filename($fasta);
    my $cmd = File::Spec->catdir(($P->{SCRIPTDIR}, "fasta_reheader.pl"));
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fasta, 1],
            ["", $name . ".filtered.fasta", 0]
        ],
        inputs => [
            $fasta,
            $validated
        ],
        outputs => [
            File::Spec->catdir(($staging_dir, $name . ".filtered.fasta"))
        ],
        description => "Filtering fasta file..."
    );
}

sub create_gff_generation_job {
    my $genome = shift;
    my $validated = shift;
    my @path = ($P->{SCRIPTDIR}, "coge_gff.pl");
    my $cmd = File::Spec->catdir(@path);
    my $org_name = sanitize_organism_name($genome->organism->name);
    my $name = "$org_name-1-name-0-0-id-" . $genome->id . "-1.gff";

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-f', $name, 0],
            ['-download_dir', '.', 0],
            ['-gid', $genome->id, 0],
            ['-upa', 1, 0],
            ['-id_type', "name", 0],
            ['-cds', 0, 0],
            ['-annos', 0, 0],
            ['-nu', 1, 0],
            ['-config', $config, 1],
        ],
        inputs => [
            $config,
            $validated
        ],
        outputs => [
            File::Spec->catdir(($CACHE, "$gid/gff", $name))
        ],
        description => "Generating gff..."
    );
}

sub create_cutadapt_job {
    my $fastq = shift;
    my $validated = shift;
    my $name = to_filename($fastq);
    my $cmd = $P->{CUTADAPT};
    die "ERROR: CUTADAPT is not in the config." unless ($cmd);

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-q', 25, 0],
            ['--quality-base=64', '', 0],
            ['-m', 17, 0],
            ['', $fastq, 1],
            ['-o', $name . '.trimmed.fastq', 1],
        ],
        inputs => [
            $fastq,
            $validated
        ],
        outputs => [
            File::Spec->catdir(($staging_dir, $name . '.trimmed.fastq'))
        ],
        description => "Running cutadapt..."
    );
}

sub create_gmap_index_job {
    my $fasta = shift;
    my $name = to_filename($fasta);
    my $cmd = $P->{GMAP_BUILD};
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
            [File::Spec->catdir(($staging_dir, $name . "-index")), 1]
        ],
        description => "Generating index..."
    );
}

sub create_gsnap_job {
    my ($fastq, $gmap) = @_;
    my $name = basename($gmap);
    my $cmd = $P->{GSNAP};
    die "ERROR: GSNAP is not in the config." unless ($cmd);

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["-D", ".", 0],
            ["-d", $name, 0],
            ["--nthreads=32", '', 0],
            ["-n", 5, 0],
            ["--format=sam", '', 0],
            ["-Q", '', 0],
            ["--gmap-mode=none", '', 0],
            ["--nofails", $fastq, 1],
            [">", $name . ".sam", 1]
        ],
        inputs => [
            $fastq,
            [$gmap, 1]
        ],
        outputs => [
            File::Spec->catdir(($staging_dir, $name . ".sam"))
        ],
        description => "Running gsnap..."
    );
}

sub create_samtools_bam_job {
    my $samfile = shift;
    my $name = to_filename($samfile);
    my $cmd = $P->{SAMTOOLS};
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
            File::Spec->catdir(($staging_dir, $name . ".bam"))
        ],
        description => "Generating bam file..."
    );
}

sub create_samtools_sort_job {
    my $bam = shift;
    my $name = to_filename($bam);
    my $cmd = $P->{SAMTOOLS};
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
            File::Spec->catdir(($staging_dir, $name . "-sorted.bam"))
        ],
        description => "Sorting bam file..."
    );
}

sub create_cufflinks_job {
    my ($gff, $fasta, $bam, undef) = @_;
    my $cmd = $P->{CUFFLINKS};
    die "ERROR: CUFFLINKS is not in the config." unless ($cmd);

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['--GTF', $gff, 1],
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
            File::Spec->catdir(($staging_dir, "genes.fpkm_tracking"))
        ],
        description => "Running cufflinks..."
    );
}

sub create_parse_cufflinks_job {
    my $cufflinks = shift;
    my $name = to_filename($cufflinks);
    my $cmd = $P->{PYTHON};
    my $script = $P->{PARSE_CUFFLINKS};
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
            File::Spec->catdir(($staging_dir, $name . ".csv"))
        ],
        description => "Parsing cufflinks to csv..."
    );
}

sub create_load_csv_job {
    my ($csv, $user) = @_;
    my $cmd = File::Spec->catdir(($P->{SCRIPTDIR}, "load_experiment.pl"));
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', qq{"$name"}, 0],
            ['-desc', qq{"Expression measurements"}, 0],
            ['-version', qq{"$version"}, 0],
            ['-restricted', $restricted, 0],
            ['-gid', $gid, 0],
            ['-source_name', qq{"$source_name"}, 0],
            ['-staging_dir', "./csv", 0],
            ['-file_type', "csv", 0],
            ['-data_file', "$csv", 0],
            ['-config', $config, 1]
        ],
        inputs => [
            $config,
            $csv
        ],
        outputs => [
            [File::Spec->catdir(($staging_dir, "csv")), 1],
            File::Spec->catdir(($staging_dir, "csv/log.done")),
        ],
        description => "Loading expression data..."
    );
}

sub create_load_bam_job {
    my ($csv, $bam, $user) = @_;
    my $cmd = File::Spec->catdir(($P->{SCRIPTDIR}, "load_experiment.pl"));
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', qq{"$name"}, 0],
            ['-desc', qq{"Mapped reads"}, 0],
            ['-version', qq{"$version"}, 0],
            ['-restricted', $restricted, 0],
            ['-gid', $gid, 0],
            ['-source_name', qq{"$source_name"}, 0],
            ['-staging_dir', "./bam", 0],
            ['-file_type', "bam", 0],
            ['-data_file', "$bam", 0],
            ['-config', $config, 1]
        ],
        inputs => [
            $config,
            $bam,
            $csv
        ],
        outputs => [
            [File::Spec->catdir(($staging_dir, "bam")), 1],
            File::Spec->catdir(($staging_dir, "csv/log.done")),
        ],
        description => "Loading mapped reads..."
    );
}

sub create_notebook_job {
    my ($bam, $csv, $user) = @_;
    my $cmd = File::Spec->catdir(($P->{SCRIPTDIR}, "create_notebook.pl"));
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-uid', $user->id, 0],
            ['-name', qq{"$name"}, 0],
            ['-desc', qq{"$description"}, 0],
            ['-page', '""', 0],
            ['-type', 2, 0],
            ['-restricted', $restricted, 0],
            ['-config', $config, 1],
            ['-log', File::Spec->catdir(($staging_dir, "log.txt")), 0],
            ['', "bam/log.txt", 0],
            ['', "csv/log.txt", 0],
        ],
        inputs => [
            $config,
            $bam,
            $csv,
            File::Spec->catdir(($staging_dir, "csv/log.done")),
            File::Spec->catdir(($staging_dir, "bam/log.done")),
        ],
        outputs => [
        ],
        description => "Creating notebook..."
    );
}

main;
