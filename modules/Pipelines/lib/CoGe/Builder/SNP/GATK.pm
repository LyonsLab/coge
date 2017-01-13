package CoGe::Builder::SNP::GATK;

use Moose;
extends 'CoGe::Builder::SNP::SNPFinder';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;


my $PICARD;
sub BUILD { # called immediately after constructor
    my $self = shift;
    $PICARD = $self->conf->{PICARD};
    unless ($PICARD) {
        CoGe::Exception::Generic->throw(message => 'Missing PICARD in config file');
    }
}

sub build {
    my $self = shift;
    my %opts = @_;
    my $bam_file = shift @{$opts{data_files}};

    my $gid = $self->request->genome->id;

    my $annotations = generate_additional_metadata();
    my @annotations2 = CoGe::Core::Metadata::to_annotations($self->params->{additional_metadata});
    push @$annotations, @annotations2;
    
    #
    # Build workflow
    #

    $self->add_task(
        $self->reheader_fasta($gid)
    );
    my $reheader_fasta = $self->previous_output;

    $self->add_task(
        $self->index_fasta($reheader_fasta)
    );

    $self->add_task(
        $self->fasta_dict($gid)
    );

    $self->add_task(
        $self->add_readgroups($bam_file)
    );

    $self->add_task_chain(
        $self->reorder_sam(
            input_fasta => $reheader_fasta,
            input_bam   => $self->previous_output
        )
    );

    $self->add_task_chain(
        $self->gatk(
            input_fasta => $reheader_fasta,
            input_bam   => $self->previous_output
        )
    );

    $self->vcf($self->previous_output);

    $self->add_task(
        $self->load_vcf(
            annotations => $annotations,
            gid         => $gid,
            vcf         => $self->vcf
        )
    );
}

sub fasta_dict {
    my $self = shift;
    my $gid = shift;

    my $fasta     = get_genome_file($gid);
    my $cache_dir = get_genome_cache_path($gid);

    my $fasta_name = to_filename($fasta);
    my $renamed_fasta = qq[$fasta_name.fa]; # mdb added 9/20/16 -- Picard expects the filename to end in .fa or .fasta
    my $fasta_dict = qq[$fasta_name.dict];
    
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

sub reorder_sam {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_bam   = $opts{input_bam};
    my $output_bam  = qq[$input_bam.reordered.bam];
    
    # mdb added 9/20/16 -- Picard expects the filename to end in .fa or .fasta
    my $fasta_name = to_filename($input_fasta);
    my $renamed_fasta = qq[$fasta_name.fa];
    
    my $done_file = qq[$output_bam.reorder.done];

    return {
        cmd => qq[ln -s $input_fasta $renamed_fasta && ln -s $input_fasta.dict $renamed_fasta.dict && java -jar $PICARD ReorderSam REFERENCE=$renamed_fasta INPUT=$input_bam OUTPUT=$output_bam CREATE_INDEX=true VERBOSITY=ERROR && touch $done_file],
        args => [],
        inputs => [
            $input_fasta,
            $input_bam
        ],
        outputs => [
            $output_bam,
            $done_file
        ],
        description => "Reorder bam file"
    };
}

sub add_readgroups {
    my $self = shift;
    my $input_bam  = shift;
    my $output_bam = qq[$input_bam.readgroups.bam];
    my $done_file  = qq[$output_bam.readgroups.done];

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

sub gatk {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_bam   = $opts{input_bam};

    my $output_vcf = 'snps.flt.vcf';

    my $params      = $self->params->{snp_params};
    my $stand_call_conf = $params->{'-stand_call_conf'} // 30;
    my $stand_emit_conf = $params->{'-stand_emit_conf'} // 10;

    my $fasta_index = qq[$input_fasta.fai];
    
    # mdb added 9/20/16 -- GATK expects the filename to end in .fa or .fasta
    my $fasta_name = to_filename($input_fasta);
    my $renamed_fasta = qq[$fasta_name.fa];
    
    my $GATK = $self->conf->{GATK};
    unless ($GATK) {
        CoGe::Exception::Generic->throw(message => 'Missing GATK in config file');
    }

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
            catfile($self->staging_dir, $output_vcf),
        ],
        description => "Identifying SNPs using GATK method"
    };
}

sub generate_additional_metadata {
    my $self = shift;
    my $params = $self->params->{snp_params};
    my $stand_call_conf = $params->{'-stand_call_conf'} || 30;
    my $stand_emit_conf = $params->{'-stand_emit_conf'} || 10; 
    
    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};
    push @annotations, qq{note|SNPs generated using GATK method, -stand_call_conf=$stand_call_conf, -stand_emit_conf=$stand_emit_conf};
    return \@annotations;
}

1;
