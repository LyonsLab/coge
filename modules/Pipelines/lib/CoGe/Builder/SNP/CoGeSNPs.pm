package CoGe::Builder::SNP::CoGeSNPs;

use v5.14;
use strict;
use warnings;

use Carp;
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
use CoGe::Accessory::Web qw(get_defaults get_job schedule_job);
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Builder::CommonTasks;

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw( build run );
}

our $CONF = CoGe::Accessory::Web::get_defaults();
our $FASTA_CACHE_DIR;

sub run {
    my %opts = @_;
    my $user = $opts{user};
    my $genome = $opts{genome};
    my $input_file = $opts{input_file};
    my $metadata = $opts{metadata};
    croak "Missing parameters" unless ($user and $genome and $input_file and $metadata);

    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Running the SNP-finder pipeline', init => 1 );

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Build the workflow
    my @jobs = build({
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        user => $user,
        wid  => $workflow->id,
        genome => $genome,
        input_file => $input_file,
        metadata => $metadata,
    });
    $workflow->add_jobs(\@jobs);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub build {
    my $opts = shift;

    # Required arguments
    my $genome = $opts->{genome};
    my $input_file = $opts->{input_file}; # path to bam file
    my $user = $opts->{user};
    my $wid = $opts->{wid};
    my $staging_dir = $opts->{staging_dir};
    my $result_dir = $opts->{result_dir};
    my $metadata = $opts->{metadata};

    my $gid = $genome->id;
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta =  to_filename($fasta_file) . ".reheader.fasta";

    $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;

    # Setup the jobs
    my @jobs;
    push @jobs, create_fasta_reheader_job(fasta => $fasta_file, reheader_fasta => $reheader_fasta, cache_dir => $FASTA_CACHE_DIR);
    push @jobs, create_fasta_index_job(fasta => $reheader_fasta, cache_dir => $FASTA_CACHE_DIR);
    push @jobs, create_samtools_job($reheader_fasta, $gid, $input_file, $staging_dir);
    push @jobs, create_load_vcf_job({
        username => $user->name,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        wid => $wid,
        gid => $gid,
        vcf => catfile($staging_dir, 'snps.vcf'),
        metadata => $metadata
    });
    
    return wantarray ? @jobs : \@jobs;
}

sub create_samtools_job {
    my ($fasta, $gid, $bam, $staging_dir) = @_;

    my $samtools = $CONF->{SAMTOOLS};
    die "ERROR: SAMTOOLS not specified in config" unless $samtools;
    my $script = catfile($CONF->{SCRIPTDIR}, 'pileup_SNPs.pl');
    die "ERROR: SCRIPTDIR not specified in config" unless $script;
    my $output_name = 'snps.vcf';

    return {
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
        description => "Identifying SNPs using the CoGe method ..."
    };
}

#sub create_load_experiment_job { # FIXME mdb 11/19/14 - replace with version in CommonTasks
#    my ($vcf, $user, $experiment, $wid) = @_;
#
#    my $cmd = catfile(($CONF->{SCRIPTDIR}, "load_experiment.pl"));
#    my $output_path = catdir($staging_dir, "load_experiment");
#
#    # Set metadata for the pipeline being used
#    my $annotations = generate_experiment_metadata();
#
#    return (
#        cmd => $cmd,
#        script => undef,
#        args => [
#            ['-user_name', $user->name, 0],
#            ['-name', '"'.$experiment->name.' (SNPs)'.'"', 0],
#            ['-desc', qq{"Single nucleotide polymorphisms"}, 0],
#            ['-version', $experiment->version, 0],
#            ['-restricted', $experiment->restricted, 0],
#            ['-gid', $experiment->genome->id, 0],
#            ['-wid', $wid, 0],
#            ['-source_name', '"'.$experiment->source->name.'"', 0],
#            ['-types', qq{"SNP"}, 0],
#            ['-annotations', $annotations, 0],
#            ['-staging_dir', "./load_experiment", 0],
#            ['-file_type', "vcf", 0],
#            ['-data_file', "$vcf", 0],
#            ['-config', $CONF->{_CONFIG_PATH}, 1]
#        ],
#        inputs => [
#            $CONF->{_CONFIG_PATH},
#            $vcf
#        ],
#        outputs => [
#            [$output_path, 1],
#            catfile($output_path, "log.done"),
#        ],
#        description => "Load SNPs as new experiment ..."
#    );
#}

sub generate_experiment_metadata {
    my @annotations = (
        qq{http://genomevolution.org/wiki/index.php/Identifying_SNPs||note|Generated by CoGe's SNP-finder Pipeline (CoGe method)},
        qq{note|Read depth generated by samtools mpileup},
        qq{note|Minimum read depth of 10},
        qq{note|Minimum high-quality (PHRED >= 20) allele count of 4},
        qq{note|Minimum allele frequency of 10%}
    );
    return join(';', @annotations);
}

1;
