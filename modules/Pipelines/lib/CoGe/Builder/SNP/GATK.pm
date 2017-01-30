package CoGe::Builder::SNP::GATK;

use v5.14;
use warnings;
use strict;

use Carp;
use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_genome_cache_path);
use CoGe::Core::Metadata qw(to_annotations);
use CoGe::Builder::CommonTasks;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(build);
our $CONF = CoGe::Accessory::Web::get_defaults();

sub build {
    my $opts = shift;

    # Required arguments
    my $genome      = $opts->{genome};
    my $input_file  = $opts->{input_file}; # path to bam file
    my $user        = $opts->{user};
    my $wid         = $opts->{wid};
    my $metadata    = $opts->{metadata};
    my $additional_metadata = $opts->{additional_metadata};
    my $params      = $opts->{params};

    # Setup paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $fasta_cache_dir = get_genome_cache_path($genome->id);
    my $fasta_file = get_genome_file($genome->id);
    my $reheader_fasta = to_filename($fasta_file) . ".reheader.faa";
    my $reheader_fasta_path = catfile($fasta_cache_dir, $reheader_fasta);

    my $annotations = generate_additional_metadata($params);
    my @annotations2 = CoGe::Core::Metadata::to_annotations($additional_metadata);
    push @$annotations, @annotations2;
    
    # Get experiment input file
    my $processed_bam_file = to_filename($input_file) . '.processed.bam';
    my $processed_bam_file2 = to_filename($input_file) . '.processed2.bam';
    my $output_vcf_file = "$processed_bam_file2.flt.vcf";

    # Build all the jobs -- TODO create cache for processed bam files
    my @tasks = (
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
            output_vcf => catfile($staging_dir, $output_vcf_file),
            params => $params
        })
    );
    my $load_vcf_task = create_load_vcf_job({
        method      => 'GATK',
        username    => $user->name,
        metadata    => $metadata,
        staging_dir => $staging_dir,
        annotations => $annotations,
        wid         => $wid,
        gid         => $genome->id,
        vcf         => catfile($staging_dir, $output_vcf_file)
    });
    push @tasks, $load_vcf_task;

    return {
        tasks => \@tasks,
        metadata => $annotations,
        done_files => [ $load_vcf_task->{outputs}->[1] ]
    };
}

sub create_fasta_dict_job {
    my $opts = shift;

    # Required arguments
    my $fasta     = $opts->{fasta};
    my $cache_dir = $opts->{cache_dir};

    my $fasta_name = to_filename($fasta);
    my $renamed_fasta = qq[$fasta_name.fa]; # mdb added 9/20/16 -- Picard expects the filename to end in .fa or .fasta
    my $fasta_dict = qq[$fasta_name.dict];
    
    my $PICARD = $CONF->{PICARD};

    return {
        cmd => qq[ln -s $fasta $renamed_fasta && java -jar $PICARD CreateSequenceDictionary REFERENCE=$renamed_fasta OUTPUT=$fasta_dict],
        args => [],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($cache_dir, $fasta_dict),
        ],
        description => "Generate fasta dictionary"
    };
}

sub create_reorder_sam_job {
    my $opts = shift;

    # Required arguments
    my $input_fasta = $opts->{input_fasta};
    my $input_bam   = $opts->{input_bam};
    my $output_bam  = $opts->{output_bam};
    
    # mdb added 9/20/16 -- Picard expects the filename to end in .fa or .fasta
    my $fasta_name = to_filename($input_fasta);
    my $renamed_fasta = qq[$fasta_name.fa];
    
    my $done_file = qq[$output_bam.reorder.done];

    my $PICARD = $CONF->{PICARD};

    return {
        cmd => qq[ln -s $input_fasta $renamed_fasta && ln -s $input_fasta.dict $renamed_fasta.dict && java -jar $PICARD ReorderSam REFERENCE=$renamed_fasta INPUT=$input_bam OUTPUT=$output_bam CREATE_INDEX=true VERBOSITY=ERROR && touch $done_file],
        args => [],
        inputs => [
            $input_fasta,
            $input_bam,
            "$input_bam.addRG.done" # from create_add_readgroups_job
        ],
        outputs => [
            $output_bam,
            $done_file
        ],
        description => "Reorder bam file"
    };
}

sub create_add_readgroups_job {
    my $opts = shift;

    # Required arguments
    my $input_bam  = $opts->{input_bam};
    my $output_bam = $opts->{output_bam};
    
    my $done_file = qq[$output_bam.addRG.done];

    my $PICARD = $CONF->{PICARD};

    return {
        cmd => qq[java -jar $PICARD AddOrReplaceReadGroups I=$input_bam O=$output_bam RGID=none RGLB=none RGPL=none RGSM=none RGPU=none && touch $done_file],
        args => [],
        inputs => [
            $input_bam
        ],
        outputs => [
            $output_bam,
            $done_file
        ],
        description => "Set read groups in bam file"
    };
}

sub create_gatk_job {
    my $opts = shift;

    # Required arguments
    my $input_fasta = $opts->{input_fasta};
    my $input_bam   = $opts->{input_bam};
    my $output_vcf  = $opts->{output_vcf};
    my $params      = $opts->{params};
    my $stand_call_conf = $params->{'-stand_call_conf'} || 30;
    my $stand_emit_conf = $params->{'-stand_emit_conf'} || 10;    

    my $fasta_index = qq[$input_fasta.fai];
    
    # mdb added 9/20/16 -- GATK expects the filename to end in .fa or .fasta
    my $fasta_name = to_filename($input_fasta);
    my $renamed_fasta = qq[$fasta_name.fa];
    
    my $GATK = $CONF->{GATK};

    return {
        cmd => qq[ln -s $input_fasta $renamed_fasta && ln -s $input_fasta.fai $renamed_fasta.fai && ln -s $input_fasta.dict $fasta_name.dict && java -jar $GATK -T HaplotypeCaller --genotyping_mode DISCOVERY --filter_reads_with_N_cigar --fix_misencoded_quality_scores],
        args =>  [
            ["-R", $renamed_fasta, 0],
            ["-I", $input_bam, 0],
            ["-stand_emit_conf", $stand_emit_conf, 0],
            ["-stand_call_conf", $stand_call_conf, 0],
            ["-o", $output_vcf, 1]
        ],
        inputs => [
            $input_bam,
            "$input_bam.reorder.done", # from create_reorder_sam_job
            $input_fasta,
            $fasta_index,
        ],
        outputs => [
            $output_vcf,
        ],
        description => "Identifying SNPs using GATK method"
    };
}

sub generate_additional_metadata {
    my $params = shift;
    my $stand_call_conf = $params->{'-stand_call_conf'} || 30;
    my $stand_emit_conf = $params->{'-stand_emit_conf'} || 10; 
    
    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};
    push @annotations, qq{note|SNPs generated using GATK method, -stand_call_conf=$stand_call_conf, -stand_emit_conf=$stand_emit_conf};
    return \@annotations;
}

1;
