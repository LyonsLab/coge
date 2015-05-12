package CoGe::Builder::SNP::Platypus;

use v5.14;
use warnings;
use strict;

use Carp;
use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Builder::CommonTasks;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(build run);
our $CONF = CoGe::Accessory::Web::get_defaults();
our $JEX = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );

sub run {
    my %opts = @_;
    my $user = $opts{user};
    my $genome = $opts{genome};
    my $input_file = $opts{input_file};
    my $metadata = $opts{metadata};
    croak "Missing parameters" unless ($user and $genome and $input_file and $metadata);

    # Create the workflow
    my $workflow = $JEX->create_workflow( name => 'Running the Playtpus SNP-finder pipeline', init => 1 );
    return unless ($workflow && $workflow->id);
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
    my $skipAnnotations = $opts->{skipAnnotations};

    # Setup paths
    my $gid = $genome->id;
    my $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta =  to_filename($fasta_file) . ".reheader.faa";

    my $annotations = generate_additional_metadata();

    my $conf = {
        staging_dir => $staging_dir,
        result_dir  => $result_dir,

        bam         => $input_file,
        fasta       => catfile($FASTA_CACHE_DIR, $reheader_fasta),
        vcf         => catfile($staging_dir, qq[snps.vcf]),

        annotations => ($skipAnnotations ? '' : join(';', @$annotations)),
        username    => $user->name,
        metadata    => $metadata,
        wid         => $wid,
        gid         => $gid
    };

    # Build the workflow's tasks
    my @tasks;
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta_file,
        reheader_fasta => $reheader_fasta,
        cache_dir => $FASTA_CACHE_DIR
    );

    push @tasks, create_fasta_index_job(
        fasta => catfile($FASTA_CACHE_DIR, $reheader_fasta),
        cache_dir => $FASTA_CACHE_DIR,
    );
    
    push @tasks, create_platypus_job($conf);
    
    my $load_vcf_task = create_load_vcf_job($conf);
    push @tasks, $load_vcf_task;

    # Save outputs for retrieval by downstream tasks
    my @done_files = (
        $load_vcf_task->{outputs}->[1]
    );
    
    my %results = (
        metadata => $annotations,
        done_files => \@done_files
    );

    return (\@tasks, \%results);
}

sub create_platypus_job {
    my $opts = shift;
#    print STDERR "create_platypus_job ", Dumper $opts, "\n";

    # Required arguments
    my $fasta = $opts->{fasta};
    my $bam = $opts->{bam};
    my $vcf = $opts->{vcf};

    my $fasta_index = qq[$fasta.fai];
    my $PLATYPUS = $CONF->{PLATYPUS} || "Platypus.py";

    return {
        cmd => qq[$PLATYPUS callVariants],
        args =>  [
            ["--bamFiles", $bam, 0],
            ["--refFile", $fasta, 0],
            ["--output", $vcf, 1],
            ["--verbosity", 0, 0],
        ],
        inputs => [
            $bam,
            $bam . '.bai',
            $fasta,
            $fasta_index,
        ],
        outputs => [
            $vcf,
        ],
        description => "Identifying SNPs using Platypus method ..."
    };
}

sub generate_additional_metadata {
    return [ qq{note|SNPs generated using Platypus method} ];
}

1;
