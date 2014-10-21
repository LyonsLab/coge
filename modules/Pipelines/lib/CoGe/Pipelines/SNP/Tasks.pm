package CoGe::Pipelines::SNP::Tasks;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(create_fasta_reheader_job create_fasta_index_job
                 create_load_vcf_job create_find_snps_job
                 create_filter_snps_job create_platypus_job);

our $CONFIG = CoGe::Accessory::Web::get_defaults();

sub create_fasta_reheader_job {
    my $opts = shift;

    # Required arguments
    my $fasta = $opts->{fasta};
    my $reheader_fasta = $opts->{reheader_fasta};
    my $cache_dir = $opts->{cache_dir};

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
    my $fasta = $opts->{fasta};
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

sub create_load_vcf_job {
    my $opts = shift;

    # Required arguments
    my $experiment = $opts->{experiment};
    my $username = $opts->{username};
    my $source_name = $opts->{source_name};
    my $staging_dir = $opts->{staging_dir};
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

sub create_platypus_job {
    my $opts = shift;

    # Required arguments
    my $reference = $opts->{fasta};
    my $alignment = $opts->{bam};
    my $vcf = $opts->{vcf};

    my $index = qq[$reference.fai];
    my $PLATYPUS = $CONFIG->{PLATYPUS} || "Platypus.py";

    return {
        cmd => qq[$PLATYPUS callVariants],
        args =>  [
            ["--bamFiles", $alignment, 0],
            ["--refFile", $reference, 0],
            ["--output", $vcf, 1],
            ["--verbosity", 0, 0],
        ],
        inputs => [
            $alignment,
            $reference,
            $index,
        ],
        outputs => [
            $vcf,
        ],
        description => "Finding SNPS using Platypus method ..."
    };
}

sub create_find_snps_job {
    my $opts = shift;

    # Required arguments
    my $reference = $opts->{fasta};
    my $alignment = $opts->{bam};
    my $snps = $opts->{bcf};

    my $subopts = {
        samtools => {
            command => $CONFIG->{SAMTOOLS} || "samtools",
            subtask => "mpileup",
            args    => {
                u => [],
                f => [],
            },
            inputs => [
                basename($reference),
                basename($alignment),
            ],
        },
        bcf => {
            command => $CONFIG->{BCFTOOLS} || "bcftools",
            subtask => "view",
            args    => {
                b => [],
                v => [],
                c => [],
                g => [],
            },
            inputs => [
            ],
        },
    };

    my @subcommands =  (
        _subcommand($subopts->{samtools}),
        _subcommand($subopts->{bcf}),
    );

    # Pipe commands together
    my $command = join " | ", @subcommands;

    # Get the output filename
    my $output = basename($snps);

    return {
        cmd => qq[$command - > $output],
        inputs => [
            $reference,
            $alignment,
        ],
        outputs => [
            $snps,
        ],
        description => "Finding SNPs using SAMtools method ...",
    };
}

sub create_filter_snps_job {
    my $opts = shift;

    # Required arguments
    my $snps = $opts->{bcf};
    my $filtered_snps = $opts->{vcf};

    # Optional arguments
    my $depth = $opts->{depth} || 100;

    my $subopts = {
        bcf => {
            command => $CONFIG->{BCFTOOLS} || "bcftools",
            subtask => "view",
            args    => {
            },
            inputs => [
                basename($snps),
            ],
        },

        vcf => {
            command => $CONFIG->{VCFTOOLS} || "vcfutils.pl",
            subtask => "varFilter",
            args    => {
                D => [$depth],
            },
            inputs => [
            ],
        },
    };

    my @subcommands =  (
        _subcommand($subopts->{bcf}),
        _subcommand($subopts->{vcf}),
    );

    # Pipe commands together
    my $command = join " | ", @subcommands;

    # Get the output filename
    my $output = basename($filtered_snps);

    return {
        cmd => qq[$command > $output],
        inputs  => [
            $snps,
        ],
        outputs => [
            $filtered_snps,
        ],
        description => "Filtering SNPs ...",
    };
}

sub _subcommand {
    my $opts = shift;
    my $cmd = $opts->{command};
    my $subtask = $opts->{subtask};
    my $args = $opts->{args};

    # create parameter string
    my @params = map { @{$args->{$_}} ? qq[-$_ @{$args->{$_}}] : qq[-$_] } keys %{$args};

    my @inputs = @{$opts->{inputs}} if $opts->{inputs};

    return qq[$cmd $subtask @params @inputs];
}

1;
