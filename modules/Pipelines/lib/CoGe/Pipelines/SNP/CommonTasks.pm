package CoGe::Pipelines::SNP::CommonTasks;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);
use Data::Dumper;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    create_fasta_reheader_job create_fasta_index_job create_load_vcf_job 
    create_bam_index_job
);

our $CONFIG = CoGe::Accessory::Web::get_defaults();

sub create_fasta_reheader_job {
    my $opts = shift;

    # Required arguments
    my $fasta          = $opts->{fasta};
    my $reheader_fasta = $opts->{reheader_fasta};
    my $cache_dir      = $opts->{cache_dir};

    my $cmd = catfile($CONFIG->{SCRIPTDIR}, "fasta_reheader.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ["", $fasta, 1],
            ["", $reheader_fasta, 0],
        ],
        inputs => [
            $fasta,
        ],
        outputs => [
            catfile($cache_dir, $reheader_fasta),
        ],
        description => "Filter fasta file...",
    };
}

sub create_fasta_index_job {
    my $opts = shift;
    
    # Required arguments
    my $fasta     = $opts->{fasta};
    my $cache_dir = $opts->{cache_dir};

    my $fasta_name = basename($fasta);
    my $fasta_index = qq[$fasta_name.fai];

    return {
        cmd => $CONFIG->{SAMTOOLS} || "samtools",
        script => undef,
        args => [
            ["faidx", $fasta, 1],
        ],
        inputs => [
            $fasta,
        ],
        outputs => [
            catfile($cache_dir, $fasta_index),
        ],
        description => "Index fasta file...",
    };
}

sub create_bam_index_job { # note: this task hasn't been tested
    my $opts = shift;
    
    # Required arguments
    my $input_bam = $opts->{input_bam};

    return {
        cmd => $CONFIG->{SAMTOOLS} || "samtools",
        script => undef,
        args => [
            ["index", $input_bam, 1],
        ],
        inputs => [
            $input_bam,
        ],
        outputs => [
            qw[$input_bam.bai]
        ],
        description => "Index bam file...",
    };
}

sub create_load_vcf_job {
    my $opts = shift;

    # Required arguments
    my $experiment = $opts->{experiment};
    my $username = $opts->{username};
    my $source_name = $opts->{source_name};
    my $staging_dir = $opts->{staging_dir};
    my $result_dir = $opts->{result_dir};
    my $annotations = $opts->{annotations};
    my $wid = $opts->{wid};
    my $gid = $opts->{gid};
    my $vcf = $opts->{vcf};

    my $cmd = catfile(($CONFIG->{SCRIPTDIR}, "load_experiment.pl"));
    my $output_path = catdir($staging_dir, "load_experiment");
    my $exp_name = $experiment->name;

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $username, 0],
            ['-name', qq["$exp_name (SNPs)"], 0],
            ['-desc', qq{"Single nucleotide polymorphisms"}, 0],
            ['-version', $experiment->version, 0],
            ['-restricted', $experiment->restricted, 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-source_name', qq[$source_name], 0],
            ['-types', qq{"SNP"}, 0],
            ['-annotations', $annotations, 0],
            ['-staging_dir', "./load_experiment", 0],
            ['-file_type', "vcf", 0],
            ['-result_dir', "'".$result_dir."'", 0],
            ['-data_file', qq[$vcf], 0],
            ['-config', $CONFIG->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            $CONFIG->{_CONFIG_PATH},
            $opts->{vcf},
        ],
        outputs => [
            [$output_path, 1],
            catfile($output_path, "log.done"),
        ],
        description => "Load SNPs as new experiment ..."
    };
}

1;
