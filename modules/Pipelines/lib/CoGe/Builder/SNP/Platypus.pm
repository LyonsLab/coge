package CoGe::Builder::SNP::Platypus;

use Moose;
extends 'CoGe::Builder::SNP::SNPFinder';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my %opts = @_;
    my $fasta_file = $opts{fasta_file};
    my ($bam_file) = @{$opts{data_files}};
    my $isSorted = $opts{is_sorted}; # input bam file is already sorted (for passing sorted bam file from Experiment.pm)

    unless ($fasta_file) {
        CoGe::Exception::Generic->throw(message => 'Missing fasta');
    }
    unless ($bam_file) {
        CoGe::Exception::Generic->throw(message => 'Missing bam');
    }

    my $gid = $self->request->genome->id;

    my $annotations = generate_additional_metadata();
    my @annotations2 = CoGe::Core::Metadata::to_annotations($self->params->{additional_metadata});
    push @$annotations, @annotations2;

    #
    # Build workflow
    #

    my $sorted_bam_file;
    if ($isSorted) {
        $sorted_bam_file = $bam_file
    }
    else {
        $self->add(
            $self->sort_bam($bam_file)
        );
        $sorted_bam_file = $self->previous_output;
        
        $self->add(
            $self->index_bam($sorted_bam_file)
        );
    }

    $self->add(
        $self->platypus(
            bam   => $sorted_bam_file,
            fasta => $fasta_file
        )
    );

    $self->vcf($self->previous_output);

    $self->add(
        $self->load_vcf(
            vcf         => $self->vcf,
            annotations => $annotations,
            gid         => $gid
        )
    );
}

sub platypus {
    my $self = shift;
    my %opts = @_;
    my $fasta = $opts{fasta};
    my $bam = $opts{bam};
    my $nCPU = 8; # number of processors to use

    my $output_file = 'snps.vcf';
    my $PLATYPUS = $self->conf->{PLATYPUS} || "Platypus.py";

    return {
        cmd => qq[export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH && $PLATYPUS callVariants], # mdb addex lib export 6/8/16
        args =>  [
            ["--bamFiles", $bam, 0],
            ["--refFile", $fasta, 0],
            ["--output", $output_file, 1],
            ["--verbosity", 0, 0],
            ["--nCPU", $nCPU, 0]
        ],
        inputs => [
            $bam,
            "$bam.bai", # bam index
            $fasta,
            "$fasta.fai" # fasta index
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Identifying SNPs using Platypus method"
    };
}

sub generate_additional_metadata {
    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};
    push @annotations, qq{note|SNPs generated using Platypus method};
    return \@annotations;
}

1;