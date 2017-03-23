package CoGe::Builder::SNP::GATK;

use Moose;
extends 'CoGe::Builder::SNP::Analyzer';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

my $JAVA_MAX_MEM = '20g';
my ($PICARD, $GATK, $NUM_CPUS);
sub BUILD { # called immediately after constructor
    my $self = shift;

    $PICARD = $self->conf->{PICARD};
    unless ($PICARD) {
        CoGe::Exception::Generic->throw(message => 'Missing PICARD in config file');
    }

    $GATK = $self->conf->{GATK};
    unless ($GATK) {
        CoGe::Exception::Generic->throw(message => 'Missing GATK in config file');
    }

    $NUM_CPUS = $self->NUM_CPUS // 16;
}

sub build {
    my $self = shift;
    my %opts = @_;
    my $fasta_file = $opts{fasta_file};
    my ($bam_file) = @{$opts{data_files}};

    unless ($fasta_file) {
        CoGe::Exception::Generic->throw(message => 'Missing fasta');
    }
    unless ($bam_file) {
        CoGe::Exception::Generic->throw(message => 'Missing bam');
    }

    my $gid = $self->request->genome->id;

    #
    # Build workflow
    #

    $self->add(
        $self->fasta_dict($fasta_file, $gid)
    );

    # Mark duplicates
    $self->add(
        $self->picard_deduplicate($bam_file)
    );

    # Add read groups -- except for BWA with -R option which sets read group values
    unless ($self->params->{alignment_params}->{tool} eq 'bwa' && $self->params->{alignment_params}->{'-R'}) {
        $self->add_to_previous(
            $self->add_readgroups($self->previous_output)
        );
    }

    # Reorder BAM
    ($bam_file) = $self->add_to_previous(
        $self->reorder_sam(
            input_fasta => $fasta_file,
            input_bam   => $self->previous_output
        )
    );

    # Realign if selected
    if ($self->params->{snp_params}->{realign}) {
        # Find intervals to analyze for realignment
        $self->add_to_previous(
            $self->gatk_RealignerTargetCreator(
                input_fasta => $fasta_file,
                input_bam   => $bam_file
            )
        );

        # Realign reads around INDELS
        ($bam_file) = $self->add_to_previous(
            $self->gatk_Realign(
                input_fasta     => $fasta_file,
                input_bam       => $bam_file,
                input_intervals => $self->previous_output
            )
        );
    }

    # Index reordered/realigned bam
    $self->add_to_previous(
        $self->index_bam($bam_file)
    );

    # Variant calling with GATK HaplotypeCaller
    $self->add_to_previous(
        $self->gatk_HaplotypeCaller(
            input_fasta => $fasta_file,
            input_bam   => $bam_file
        )
    );

    # Set pipeline output
    $self->vcf($self->previous_output);

    # Load VCF experiment
    $self->add(
        $self->load_vcf(
            annotations => $self->generate_additional_metadata(), #TODO use metadata file instead
            gid         => $gid,
            vcf         => $self->vcf
        )
    );
}

sub reorder_sam {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_bam   = $opts{input_bam};
    my $output_bam  = qq[$input_bam.reordered.bam];
    
    # mdb added 9/20/16 -- Picard expects the filename to end in .fa or .fasta
    #my $fasta_name = to_filename($input_fasta);
    #my $renamed_fasta = qq[$fasta_name.fa];
    my $renamed_fasta = qq[$input_fasta.fa];
    
    my $done_file = qq[$output_bam.reorder.done];

    return {
        cmd => qq[ln -sf $input_fasta $renamed_fasta && ln -sf $input_fasta.dict $renamed_fasta.dict && java -jar $PICARD ReorderSam REFERENCE=$renamed_fasta INPUT=$input_bam OUTPUT=$output_bam CREATE_INDEX=true VERBOSITY=ERROR && touch $done_file],
        args => [],
        inputs => [
            $input_fasta,
            $input_bam,
            qq[$input_fasta.dict]
        ],
        outputs => [
            $output_bam,
            $done_file
        ],
        description => "Reorder BAM file"
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
        description => "Set read groups in BAM file"
    };
}

sub gatk_RealignerTargetCreator {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_bam   = $opts{input_bam};

    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;

    my $snp_params = $self->params->{snp_params};
    my $filter_reads_with_N_cigar     = $snp_params->{'--filter_reads_with_N_cigar'}     ? '--filter_reads_with_N_cigar'     : '';
    my $fix_misencoded_quality_scores = $snp_params->{'--fix_misencoded_quality_scores'} ? '--fix_misencoded_quality_scores' : '';

    my $output_file = to_filename_base($input_bam) . '.realignment.intervals';

    my $fasta_index = qq[$input_fasta.fai];

    # mdb added 9/20/16 -- GATK expects the filename to end in .fa or .fasta
    my $fasta_name = to_filename($input_fasta);
    my $renamed_fasta = qq[$fasta_name.fa];

    my $args = [
        #['-U', 'ALLOW_N_CIGAR_READS', 0], # fixes "ERROR MESSAGE: Unsupported CIGAR operator N in read" for test file test_rna_seq_data_0.17M_reads.fastq.gz
        ['-R', $renamed_fasta,        0],
        ['-I', $input_bam,            0],
        ['-o', $output_file,          1]
    ];
    #push @$args, ['-fixMisencodedQuals', '', 0] if ($encoding == 64);

    return {
        cmd => "ln -sf $input_fasta $renamed_fasta && " . # GATK expects the filename to end in .fa or .fasta
               "ln -sf $input_fasta.fai $renamed_fasta.fai && " .
               "ln -sf $input_fasta.dict $fasta_name.dict && " .
               qq[java -Xmx$JAVA_MAX_MEM -jar $GATK -T RealignerTargetCreator $filter_reads_with_N_cigar $fix_misencoded_quality_scores],
        args =>  $args,
        inputs => [
            $input_bam,
            $input_fasta,
            $fasta_index,
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Finding intervals to analyze"
    };
}

sub gatk_Realign {
    my $self = shift;
    my %opts = @_;
    my $input_fasta     = $opts{input_fasta};
    my $input_bam       = $opts{input_bam};
    my $input_intervals = $opts{input_intervals};

    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;

    my $snp_params = $self->params->{snp_params};
    my $filter_reads_with_N_cigar     = $snp_params->{'--filter_reads_with_N_cigar'}     ? '--filter_reads_with_N_cigar'     : '';
    my $fix_misencoded_quality_scores = $snp_params->{'--fix_misencoded_quality_scores'} ? '--fix_misencoded_quality_scores' : '';

    my $output_file = to_filename_base($input_bam) . '.realigned.bam';

    my $fasta_index = qq[$input_fasta.fai];

    # mdb added 9/20/16 -- GATK expects the filename to end in .fa or .fasta
    my $fasta_name = to_filename($input_fasta);
    my $renamed_fasta = qq[$fasta_name.fa];

    my $args = [
        ["-R",               $renamed_fasta,   0],
        ["-I",               $input_bam,       0],
        ["-targetIntervals", $input_intervals, 0],
        ["-o",               $output_file,     1]
    ];

    return {
        cmd => "ln -sf $input_fasta $renamed_fasta && " . # GATK expects the filename to end in .fa or .fasta
               "ln -sf $input_fasta.fai $renamed_fasta.fai && " .
               "ln -sf $input_fasta.dict $fasta_name.dict && " .
                qq[java -Xmx$JAVA_MAX_MEM -jar $GATK -T IndelRealigner $filter_reads_with_N_cigar $fix_misencoded_quality_scores],
        args =>  $args,
        inputs => [
            $input_bam,
            $input_fasta,
            $fasta_index,
            $input_intervals
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Realigning reads"
    };
}

sub gatk_HaplotypeCaller {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_bam   = $opts{input_bam};

    my $params = $self->params->{snp_params};
    my $stand_call_conf   = $params->{'-stand_call_conf'};
    my $emitRefConfidence = $params->{'--emitRefConfidence'};
    my $pcr_indel_model   = $params->{'--pcr_indel_model'};
    my $filter_reads_with_N_cigar     = $params->{'--filter_reads_with_N_cigar'}     ? '--filter_reads_with_N_cigar'     : '';
    my $fix_misencoded_quality_scores = $params->{'--fix_misencoded_quality_scores'} ? '--fix_misencoded_quality_scores' : '';

    my $output_vcf = to_filename_base($input_bam) . ($emitRefConfidence ? '.g.vcf' : '.vcf');

    my $args = [
        ['-R', qq[$input_fasta.fa], 0],
        ['-I', $input_bam,          1],
        ['-o', $output_vcf,         1]
    ];

    push @$args, ['--emitRefConfidence', $emitRefConfidence,  0] if ($emitRefConfidence);
    push @$args, ['--pcr_indel_model',   $pcr_indel_model,    0] if ($pcr_indel_model);
    push @$args, ['-stand_call_conf',    $stand_call_conf,    0] if ($stand_call_conf);

    return {
        cmd => "ln -sf $input_fasta $input_fasta.fa && " . # GATK expects the filename to end in .fa or .fasta
               "ln -sf $input_fasta.fai $input_fasta.fa.fai && " .
               "ln -sf $input_fasta.dict $input_fasta.fa.dict && " .
                qq[java -Xmx$JAVA_MAX_MEM -jar $GATK -T HaplotypeCaller -nct $NUM_CPUS $filter_reads_with_N_cigar $fix_misencoded_quality_scores],
        args => $args,
        inputs => [
            $input_bam,
            qq[$input_bam.bai],
            $input_fasta,
            qq[$input_fasta.fai],
            qq[$input_fasta.dict]
        ],
        outputs => [
            catfile($self->staging_dir, $output_vcf)
        ],
        description => "Identifying SNPs using GATK HaplotypeCaller"
    };
}

sub generate_additional_metadata {
    my $self = shift;
    my $params = $self->params->{snp_params};
    my $filter_reads_with_N_cigar     = $params->{'--filter_reads_with_N_cigar'}     ? '--filter_reads_with_N_cigar'     : '';
    my $fix_misencoded_quality_scores = $params->{'--fix_misencoded_quality_scores'} ? '--fix_misencoded_quality_scores' : '';
    my $realign = $params->{'realign'} ? 'yes' : 'no';

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};
    push @annotations, qq{note|SNPs generated using GATK method realign=$realign $filter_reads_with_N_cigar $fix_misencoded_quality_scores } .
        join(' ', map { ($params->{$_} ? $_.' '.$params->{$_} : '') } ('--emitRefConfidence', '--pcr_indel_model', '-stand_call_conf'));
    return \@annotations;
}

1;
