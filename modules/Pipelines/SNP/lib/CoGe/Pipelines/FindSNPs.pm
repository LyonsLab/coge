package CoGe::Pipelines::FindSNPs;

use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Basename qw(fileparse basename dirname);
use File::Path qw(mkpath);
use File::Spec::Functions qw(catfile catdir);
use Getopt::Long qw(GetOptions);
use JSON qw(decode_json);
use URI::Escape::JavaScript qw(unescape);

use CoGe::Accessory::TDS qw(read);
use CoGe::Accessory::Workflow;
use CoGe::Accessory::Jex;
use CoGe::Core::Storage qw(get_genome_file get_experiment_files get_workflow_paths);
use CoGe::Accessory::Web qw(get_defaults get_job schedule_job);

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw( run );
}

our ($CONF, $staging_dir, $result_dir, $log_file, $METADATA, $FASTA_CACHE_DIR );

sub run {
    my %opts = @_;
    my $db = $opts{db};
    my $experiment = $opts{experiment};
    my $user = $opts{user};

    my $eid = $experiment->id;
    my $genome = $experiment->genome;

    $CONF = CoGe::Accessory::Web::get_defaults();

    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Running the SNP-finder pipeline', init => 1 );

    # Setup log file, staging, and results paths
    ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($staging_dir, 'log_main.txt') );

    $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $genome->id, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;

    my $fasta_file = get_genome_file($genome->id);
    my $files = get_experiment_files($eid, $experiment->data_type);
    my $bam_file = shift @$files;

    # Setup the jobs
    my $filtered_file = to_filename($fasta_file) . ".filtered.fasta";
    $workflow->add_job(
        create_fasta_reheader_job($fasta_file, $genome->id, $filtered_file)
    );
    $workflow->add_job(
        create_fasta_index_job($filtered_file, $genome->id)
    );
    $workflow->add_job(
        create_samtools_job($filtered_file, $genome->id, $bam_file)
    );
    $workflow->add_job(
        create_load_experiment_job(catfile($staging_dir, 'snps.vcf'), $user, $experiment)
    );

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub to_filename { # FIXME: move into Utils module
    my ($name, undef, undef) = fileparse(shift, qr/\.[^.]*/);
    return $name;
}

sub create_fasta_reheader_job {
    my ($fasta, $gid, $output) = @_;

    my $cmd = catfile($CONF->{SCRIPTDIR}, "fasta_reheader.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fasta, 1],
            ["", $output, 0]
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($FASTA_CACHE_DIR, $output)
        ],
        description => "Filter fasta file..."
    );
}

sub create_fasta_index_job {
    my ($fasta, $gid) = @_;

    my $samtools = $CONF->{SAMTOOLS};
    die "ERROR: SAMTOOLS not specified in config" unless $samtools;

    return (
        cmd => $samtools,
        script => undef,
        args => [
            ['faidx', '', 0],
            ['', $fasta, 1]
        ],
        inputs => [
            catfile($FASTA_CACHE_DIR, $fasta)
        ],
        outputs => [
            catfile($FASTA_CACHE_DIR, $fasta) . '.fai'
        ],
        description => "Index fasta file..."
    );
}

sub create_samtools_job {
    my ($fasta, $gid, $bam) = @_;

    my $samtools = $CONF->{SAMTOOLS};
    die "ERROR: SAMTOOLS not specified in config" unless $samtools;
    my $script = catfile($CONF->{SCRIPTDIR}, 'pileup_SNPs.pl');
    die "ERROR: SCRIPTDIR not specified in config" unless $script;
    my $output_name = 'snps.vcf';

    return (
        cmd => $samtools,
        script => undef,
        args => [
            ['mpileup', '', 0],
            ['-f', '', 0],
            ['', $fasta, 1],
            ['', $bam, 1],
            ['|', $script, 0],
            ['>', $output_name,  0]
        ],
        inputs => [
            catfile($FASTA_CACHE_DIR, $fasta),
            catfile($FASTA_CACHE_DIR, $fasta) . '.fai',
            $bam
        ],
        outputs => [
            catfile($staging_dir, $output_name)
        ],
        description => "Identify SNPs ..."
    );
}

sub create_load_experiment_job {
    my ($vcf, $user, $experiment) = @_;

    my $cmd = catfile(($CONF->{SCRIPTDIR}, "load_experiment.pl"));
    my $output_path = catdir($staging_dir, "load_experiment");

    # Set metadata for the pipeline being used
    my $annotations = generate_experiment_metadata();

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', '"'.$experiment->name.' (SNPs)'.'"', 0],
            ['-desc', qq{"Single nucleotide polymorphisms"}, 0],
            ['-version', $experiment->version, 0],
            ['-restricted', $experiment->restricted, 0],
            ['-gid', $experiment->genome->id, 0],
            ['-source_name', '"'.$experiment->source->name.'"', 0],
            ['-types', qq{"SNP"}, 0],
            ['-annotations', $annotations, 0],
            ['-staging_dir', "./load_experiment", 0],
            ['-result_dir', $result_dir, 0],
            ['-file_type', "vcf", 0],
            ['-data_file', "$vcf", 0],
            ['-config', $CONF->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONF->{_CONFIG_PATH},
            $vcf
        ],
        outputs => [
            [$output_path, 1],
            catfile($output_path, "log.done"),
        ],
        description => "Load SNPs as new experiment ..."
    );
}

sub generate_experiment_metadata {
    my @annotations = (
        qq{http://genomevolution.org/wiki/index.php/Identifying_SNPs||note|Generated by CoGe's SNP-finder Pipeline},
        qq{note|Read depth generated by samtools mpileup},
        qq{note|Minimum read depth of 10},
        qq{note|Minimum high-quality (PHRED >= 20) allele count of 4},
        qq{note|Minimum allele frequency of 10%}
    );
    return '"' . join(';', @annotations) . '"';
}

1;
