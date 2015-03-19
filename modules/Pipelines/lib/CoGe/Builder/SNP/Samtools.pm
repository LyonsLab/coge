package CoGe::Builder::SNP::Samtools;

use v5.14;
use warnings;
use strict;

use Carp;
use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);

use CoGe::Accessory::Jex;
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Builder::CommonTasks;

our $CONF = CoGe::Accessory::Web::get_defaults();
our $JEX = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
    require Exporter;

    $VERSION   = 0.1;
    @ISA       = qw(Exporter);
    @EXPORT    = qw(run build);
    #@EXPORT_OK = qw(build);
}

sub run {
    my %opts = @_;
    my $user = $opts{user};
    my $genome = $opts{genome};
    my $input_file = $opts{input_file};
    my $metadata = $opts{metadata};
    croak "Missing parameters" unless ($user and $genome and $input_file and $metadata);

    # Create the workflow
    my $workflow = $JEX->create_workflow( name => 'Running the SAMtools SNP-finder pipeline', init => 1 );
    my ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Build the workflow
    my @tasks = build({
        user => $user,
        wid  => $workflow->id,
        genome => $genome,
        input_file => $input_file,
        metadata => $metadata,
    });
    $workflow->add_jobs(\@tasks);

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
    my $genome = $opts->{genome};
    my $input_file = $opts->{input_file}; # path to bam file
    my $user = $opts->{user};
    my $wid = $opts->{wid};
    my $metadata = $opts->{metadata};
    my $params = $opts->{params};

    # Setup paths
    my $gid = $genome->id;
    my $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta =  to_filename($fasta_file) . ".reheader.faa";

    my $conf = {
        staging_dir => $staging_dir,
        result_dir  => $result_dir,

        bam         => $input_file,
        fasta       => catfile($FASTA_CACHE_DIR, $reheader_fasta),
        bcf         => catfile($staging_dir, qq[snps.raw.bcf]),
        vcf         => catfile($staging_dir, qq[snps.flt.vcf]),

        username    => $user->name,
        metadata    => $metadata,
        wid         => $wid,
        gid         => $gid,
        
        params      => $params
    };

    # Build the workflow's tasks
    my @tasks;
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta_file,
        reheader_fasta => $reheader_fasta,
        cache_dir => $FASTA_CACHE_DIR
    );

    push @tasks, create_find_snps_job($conf);
    
    push @tasks, create_filter_snps_job($conf);
    
    my $load_vcf_task = create_load_vcf_job($conf);
    push @tasks, $load_vcf_task;

    # Save outputs for retrieval by downstream tasks
    my @done_files = (
        $load_vcf_task->{outputs}->[1]
    );
    
    my %results = (
        metadata => generate_additional_metadata($params),
        done_files => \@done_files
    );

    return (\@tasks, \%results);
}

sub create_find_snps_job {
    my $opts = shift;

    # Required arguments
    my $reference = $opts->{fasta};
    my $alignment = $opts->{bam};
    my $snps = $opts->{bcf};

    my $subopts = {
        samtools => {
            command => $CONF->{SAMTOOLS} || "samtools",
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
            command => $CONF->{BCFTOOLS} || "bcftools",
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
        description => "Identifying SNPs using SAMtools method ...",
    };
}

sub create_filter_snps_job {
    my $opts = shift;

    # Required arguments
    my $snps = $opts->{bcf};
    my $filtered_snps = $opts->{vcf};

    # Optional arguments
    my $params = $opts->{params};
    my $min_read_depth = $params->{'min-read-depth'} || 6;
    my $max_read_depth = $params->{'max-read-depth'} || 10;

    my $subopts = {
        bcf => {
            command => $CONF->{BCFTOOLS} || "bcftools",
            subtask => "view",
            args    => {},
            inputs => [ basename($snps) ]
        },

        vcf => {
            command => $CONF->{VCFTOOLS} || "vcfutils.pl",
            subtask => "varFilter",
            args    => {
                d => [$min_read_depth],
                D => [$max_read_depth]
            },
            inputs => []
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

sub generate_additional_metadata {
    my $params = shift;
    my $min_read_depth = $params->{'min-read-depth'} || 6;
    my $max_read_depth = $params->{'max-read-depth'} || 10;
    return [ qq{note|SNPs generated using SAMtools method, min read depth $min_read_depth, max read depth $max_read_depth} ];
}

1;
