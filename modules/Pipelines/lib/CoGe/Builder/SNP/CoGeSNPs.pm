package CoGe::Builder::SNP::CoGeSNPs;

use v5.14;
use strict;
use warnings;

use Carp;
use Data::Dumper qw(Dumper);
use File::Basename qw(fileparse basename dirname);
use File::Spec::Functions qw(catfile catdir);
use Getopt::Long qw(GetOptions);
use JSON qw(decode_json);
use URI::Escape::JavaScript qw(unescape);

use CoGe::Accessory::Web qw(get_defaults get_command_path get_job schedule_job);
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_genome_cache_path);
use CoGe::Core::Metadata qw(to_annotations);
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

sub build {
    my $opts = shift;

    # Required arguments
    my $genome = $opts->{genome};
    my $input_file = $opts->{input_file}; # path to bam file
    my $user = $opts->{user};
    my $wid = $opts->{wid};
    my $metadata = $opts->{metadata};
    my $additional_metadata = $opts->{additional_metadata};
    my $params = $opts->{params};

    # Setup paths
    my $gid = $genome->id;
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta =  to_filename($fasta_file) . ".reheader.faa";
    
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);

    $FASTA_CACHE_DIR = get_genome_cache_path($gid);
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;

    # Build the workflow's tasks
    my @tasks;
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta_file, 
        reheader_fasta => $reheader_fasta, 
        cache_dir => $FASTA_CACHE_DIR
    );
    
    push @tasks, create_fasta_index_job(
        fasta => catfile($FASTA_CACHE_DIR, $reheader_fasta), 
        cache_dir => $FASTA_CACHE_DIR
    );
    
    push @tasks, create_samtools_job(
        reheader_fasta => $reheader_fasta, 
        gid => $gid, 
        input_file => $input_file, 
        staging_dir => $staging_dir
    );
    
    my $annotations = generate_additional_metadata($params);
    my @annotations2 = CoGe::Core::Metadata::to_annotations($additional_metadata);
    push @$annotations, @annotations2;
    
    my $load_vcf_task = create_load_vcf_job({
        method => 'CoGe',
        username => $user->name,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        wid => $wid,
        gid => $gid,
        vcf => catfile($staging_dir, "$input_file.vcf"),
        metadata => $metadata,
        annotations => $annotations
    });
    push @tasks, $load_vcf_task;
    
    return {
        tasks => \@tasks,
        metadata => $annotations,
        done_files => [ $load_vcf_task->{outputs}->[1] ]
    };
}

sub create_samtools_job {
    my %opts = @_;

    # Required arguments
    my $reheader_fasta = $opts{reheader_fasta};
    my $gid = $opts{gid};
    my $bam_file = $opts{input_file};
    my $staging_dir = $opts{staging_dir};
    
    # Optional arguments
    my $params = $opts{params};
    my $min_read_depth   = $params->{'min-read-depth'} || 10;
    my $min_base_quality = $params->{'min-base-quality'} || 20;
    my $min_allele_freq  = $params->{'min-allele-freq'} || 0.1;
    my $min_allele_count = $params->{'min-allele-count'} || 4;
    my $scale            = $params->{scale} || 32;
    
    my $samtools = get_command_path('SAMTOOLS');
    
    die "ERROR: SCRIPTDIR not specified in config" unless $CONF->{SCRIPTDIR};
    my $filter_script = catfile($CONF->{SCRIPTDIR}, 'pileup_SNPs.pl');
    $filter_script .= ' min_read_depth=' . $min_read_depth;
    $filter_script .= ' min_base_quality=' . $min_base_quality;
    $filter_script .= ' min_allele_freq=' . $min_allele_freq;
    $filter_script .= ' min_allele_count=' . $min_allele_count;
    $filter_script .= ' quality_scale=' . $scale;
    
    my $output_name = "$bam_file.vcf";

    return {
        cmd => $samtools,
        script => undef,
        args => [
            ['mpileup', '', 0],
            ['-f', '', 0],
            ['', $reheader_fasta, 1],
            ['', $bam_file, 1],
            ['|', $filter_script, 0],
            ['>', $output_name,  0]
        ],
        inputs => [
            catfile($FASTA_CACHE_DIR, $reheader_fasta),
            catfile($FASTA_CACHE_DIR, $reheader_fasta) . '.fai',
            $bam_file
        ],
        outputs => [
            catfile($staging_dir, $output_name)
        ],
        description => "Identifying SNPs using the CoGe method"
    };
}

sub generate_additional_metadata {
    my $params = shift;
    my $min_read_depth   = $params->{'min-read-depth'} || 10;
    my $min_base_quality = $params->{'min-base-quality'} || 20;
    my $min_allele_freq  = $params->{'min-allele-freq'} || 0.1;
    my $min_allele_count = $params->{'min-allele-count'} || 4;
    my $scale            = $params->{scale} || 32;

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};
    push @annotations, (
        qq{note|SNPs generated using CoGe method},
        qq{note|Minimum read depth of $min_read_depth},
        qq{note|Minimum high-quality (PHRED >= $min_base_quality) allele count of $min_allele_count, FASTQ encoding $scale},
        qq{note|Minimum allele frequency of } . $min_allele_freq * 100 . '%'
    );
    return \@annotations;
}

1;
