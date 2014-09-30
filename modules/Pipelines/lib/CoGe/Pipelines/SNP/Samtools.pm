package CoGe::Pipelines::SNP::Samtools;

use v5.14;
use warnings;
use strict;

use Carp;
use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(fileparse basename);

use CoGe::Core::Storage qw(get_genome_file get_experiment_files get_workflow_paths);
use CoGe::Accessory::Web;
use CoGe::Accessory::Jex;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(build run);
our $CONFIG = CoGe::Accessory::Web::get_defaults();
our $JEX = CoGe::Accessory::Jex->new( host => $CONFIG->{JOBSERVER}, port => $CONFIG->{JOBPORT} );

sub run {
    my %opts = @_;

    # Required arguments
    my $experiment = $opts{experiment} or croak "An experiment must be specified";
    my $user = $opts{user} or croak "A user was not specified";

    my $workflow = $JEX->create_workflow( name => 'Running the SNP-finder pipeline', init => 1 );
    my ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($staging_dir, 'debug.log') );

    my @jobs = build({
        experiment => $experiment,
        staging_dir => $staging_dir,
        user => $user,
        wid  => $workflow->id,
    });

    # Add all the jobs to the workflow
    foreach (@jobs) {
        $workflow->add_job(%{$_});
    }

    # Submit the workflow
    my $result = $JEX->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub build {
    my $opts = shift;

    # Required arguments
    my $experiment = $opts->{experiment};
    my $user = $opts->{user};
    my $wid = $opts->{wid};
    my $staging_dir = $opts->{staging_dir};

    my $genome = $experiment->genome;
    my $fasta_cache_dir = catdir($CONFIG->{CACHEDIR}, $genome->id, "fasta");

    my $fasta_file = get_genome_file($genome->id);
    my $files = get_experiment_files($experiment->id, $experiment->data_type);
    my $bam_file = shift @$files;
    my $basename = to_filename($bam_file);
    my $reheader_fasta =  to_filename($fasta_file) . ".filtered.fasta";

    my $conf = {
        staging_dir    => $staging_dir,

        bam            => $bam_file,
        fasta          => catfile($fasta_cache_dir, $reheader_fasta),
        bcf            => catfile($staging_dir, qq[snps.raw.bcf]),
        vcf            => catfile($staging_dir, qq[snps.flt.vcf]),

        experiment     => $experiment,
        username       => $user->name,
        source_name    => $experiment->source->name,
        wid            => $wid,
        gid            => $genome->id,
    };

    my @jobs;

    # Build all the job
    push @jobs, create_fasta_reheader_job({
        fasta => $fasta_file,
        cache_dir => $fasta_cache_dir,
        reheader_fasta => $reheader_fasta,
    });

    push @jobs, create_find_snps_job($conf);
    push @jobs, create_filter_snps_job($conf);
    push @jobs, create_load_vcf_job($conf);

    return @jobs;
}

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

sub create_load_vcf_job {
    my $opts = shift;

    # Required arguments
    my $experiment = $opts->{experiment};
    my $username = $opts->{username};
    my $source_name = $opts->{source_name};
    my $staging_dir = $opts->{staging_dir};
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
            ['-annotations', $opts->{annotations}, 0],
            ['-staging_dir', "./load_experiment", 0],
            ['-file_type', "vcf", 0],
            ['-data_file', qq[$opts->{vcf}], 0],
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

sub to_filename { # FIXME: move into Utils module
    my ($name, undef, undef) = fileparse(shift, qr/\.[^.]*/);
    return $name;
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
