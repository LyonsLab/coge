package CoGe::Builder::SNP::GATK;

use v5.14;
use warnings;
use strict;

use Carp;
use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);

use CoGe::Accessory::Jex;
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Accessory::Web;
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Builder::CommonTasks;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(build run);
our $CONFIG = CoGe::Accessory::Web::get_defaults();
our $JEX = CoGe::Accessory::Jex->new( host => $CONFIG->{JOBSERVER}, port => $CONFIG->{JOBPORT} );

sub run {
    my %opts = @_;
    my $user = $opts{user};
    my $genome = $opts{genome};
    my $input_file = $opts{input_file};
    my $metadata = $opts{metadata};
    croak "Missing parameters" unless ($user and $genome and $input_file and $metadata);

    # Create a workflow
    my $workflow = $JEX->create_workflow( name => 'Running the GATK SNP-finder pipeline', init => 1 );
    return unless ($workflow && $workflow->id);
    my ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Build the workflow
    my @jobs = build({
        user => $user,
        wid  => $workflow->id,
        genome => $genome,
        input_file => $input_file,
        metadata => $metadata,
    });
    $workflow->add_jobs(\@jobs);

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
    my $user        = $opts->{user};
    my $wid         = $opts->{wid};
    my $metadata    = $opts->{metadata};

    # Get genome and associated files/paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $fasta_cache_dir = catdir($CONFIG->{CACHEDIR}, $genome->id, "fasta");
    my $fasta_file = get_genome_file($genome->id);
    my $reheader_fasta = to_filename($fasta_file) . ".reheader.faa";
    my $reheader_fasta_path = catfile($fasta_cache_dir, $reheader_fasta);
    
    # Get experiment input file
    my $processed_bam_file = to_filename($input_file) . '.processed.bam';
    my $processed_bam_file2 = to_filename($input_file) . '.processed2.bam';
    my $output_vcf_file = qq[snps.flt.vcf];

    # Build all the jobs -- TODO create cache for processed bam files
    my @jobs = (
        create_fasta_reheader_job(
            fasta => $fasta_file,
            reheader_fasta => $reheader_fasta,
            cache_dir => $fasta_cache_dir,
        ),
        create_fasta_index_job(
            fasta => $reheader_fasta_path,
            cache_dir => $fasta_cache_dir
        ),
        create_fasta_dict_job({
            fasta => $reheader_fasta_path,
            cache_dir => $fasta_cache_dir
        }),
        create_add_readgroups_job({
            input_bam => $input_file,
            output_bam => catfile($staging_dir, $processed_bam_file),
        }),
        create_reorder_sam_job({
            input_fasta => $reheader_fasta_path,
            input_bam => catfile($staging_dir, $processed_bam_file),
            output_bam => catfile($staging_dir, $processed_bam_file2),            
        }),
        create_gatk_job({
            input_fasta => $reheader_fasta_path,
            input_bam => catfile($staging_dir, $processed_bam_file2),
            output_vcf => catfile($staging_dir, $output_vcf_file)
        }),
        create_load_vcf_job({
            method      => 'GATK',
            username    => $user->name,
            metadata    => $metadata,
            staging_dir => $staging_dir,
            result_dir  => $result_dir,
            annotations => generate_experiment_metadata(),
            wid         => $wid,
            gid         => $genome->id,
            vcf         => $output_vcf_file
        })
    );

    return wantarray ? @jobs : \@jobs;
}

sub create_fasta_dict_job {
    my $opts = shift;

    # Required arguments
    my $fasta     = $opts->{fasta};
    my $cache_dir = $opts->{cache_dir};

    my $fasta_name = to_filename($fasta);
    my $fasta_dict = qq[$fasta_name.dict];
    
    my $PICARD = $CONFIG->{PICARD};

    return {
        cmd => qq[java -jar $PICARD CreateSequenceDictionary REFERENCE=$fasta OUTPUT=$fasta_dict],
        args => [],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($cache_dir, $fasta_dict),
        ],
        description => "Generate fasta dictionary ..."
    };
}

sub create_reorder_sam_job {
    my $opts = shift;

    # Required arguments
    my $input_fasta = $opts->{input_fasta};
    my $input_bam   = $opts->{input_bam};
    my $output_bam  = $opts->{output_bam};

    my $PICARD = $CONFIG->{PICARD};

    return {
        cmd => qq[java -jar $PICARD ReorderSam REFERENCE=$input_fasta INPUT=$input_bam OUTPUT=$output_bam CREATE_INDEX=true VERBOSITY=ERROR],
        args => [],
        inputs => [
            $input_fasta,
            $input_bam
        ],
        outputs => [
            $output_bam
        ],
        description => "Reorder bam file ..."
    };
}

sub create_add_readgroups_job {
    my $opts = shift;

    # Required arguments
    my $input_bam  = $opts->{input_bam};
    my $output_bam = $opts->{output_bam};

    my $PICARD = $CONFIG->{PICARD};

    return {
        cmd => qq[java -jar $PICARD AddOrReplaceReadGroups I=$input_bam O=$output_bam RGID=none RGLB=none RGPL=none RGSM=none RGPU=none],
        args => [],
        inputs => [
            $input_bam
        ],
        outputs => [
            $output_bam,
        ],
        description => "Set read groups in bam file ..."
    };
}

sub create_gatk_job {
    my $opts = shift;

    # Required arguments
    my $input_fasta = $opts->{input_fasta};
    my $input_bam   = $opts->{input_bam};
    my $output_vcf  = $opts->{output_vcf};

    my $fasta_index = qq[$input_fasta.fai];
    my $GATK = $CONFIG->{GATK};

    return {
        cmd => qq[java -jar $GATK -T HaplotypeCaller --genotyping_mode DISCOVERY --filter_reads_with_N_cigar --fix_misencoded_quality_scores],
        args =>  [
            ["-R", $input_fasta, 0],
            ["-I", $input_bam, 0],
            ["-stand_emit_conf", 10, 0],
            ["-stand_call_conf", 30, 0],
            ["-o", $output_vcf, 0]
        ],
        inputs => [
            $input_bam,
            $input_fasta,
            $fasta_index,
        ],
        outputs => [
            $output_vcf,
        ],
        description => "Identifying SNPs using GATK method ..."
    };
}

sub generate_experiment_metadata {
    my @annotations = (
        qq{http://genomevolution.org/wiki/index.php/Identifying_SNPs||note|Generated by CoGe's SNP-finder Pipeline (GATK method)},
    );
    return '"' . join(';', @annotations) . '"';
}

1;
