package CoGe::Builder::Alignment::Bismark;

use Moose;
extends 'CoGe::Builder::Alignment::Aligner';

use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage qw(get_genome_file get_genome_cache_path);
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my %opts = @_;
    my $fasta = $opts{fasta_file}; # reheadered fasta file
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    # Add index task
    $self->add(
        $self->bismark_index($fasta)
    );

    # Add alignment task
    $self->add_to_previous(
        $self->bismark_alignment($fastq)
    );
    push @{$self->bam}, $self->previous_output;
}

sub bismark_index {
    my $self = shift;
    my $fasta = shift;

    my $gid = $self->request->genome->id;
    my $cache_dir = catdir(get_genome_cache_path($gid), "bismark_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');

    my $done_file = 'bismark_genome_preparation.done';

    my $cmd = $self->conf->{BISMARK_DIR} ? catfile($self->conf->{BISMARK_DIR}, 'bismark_genome_preparation') : 'bismark_genome_preparation';
    $cmd = "cp $fasta $name.fa && " . # bismark requires fasta file to end in .fa or .fasta, not .faa
           "nice $cmd $cache_dir && " .
           "touch $done_file";

    $self->index([
        catfile($cache_dir, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.1.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.2.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.3.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.4.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.rev.1.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'CT_conversion', 'BS_CT.rev.2.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.1.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.2.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.3.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.4.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.rev.1.bt2'),
        catfile($cache_dir, 'Bisulfite_Genome', 'GA_conversion', 'BS_GA.rev.2.bt2'),
    ]);

    return {
        cmd         => $cmd,
        args        => [],
        inputs      => [
            $fasta
        ],
        outputs     => [
            @{$self->index},
            catfile($cache_dir, $done_file)
        ],
        description => "Indexing genome sequence with Bismark"
    };
}

sub bismark_alignment {
    my $self  = shift;
    my $fastq = shift;

    my $gid         = $self->request->genome->id;
    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;
    my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{alignment_params} // {};
    my $N = $alignment_params->{'-N'} // 0;
    my $L = $alignment_params->{'-L'} // 20;

    my $index_path = catdir(get_genome_cache_path($gid), "bismark_index");

    # Build up command/arguments string
    my $cmd = $self->conf->{BISMARK_DIR} ? catfile($self->conf->{BISMARK_DIR}, 'bismark') : 'bismark';
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ['-p', 4, 0], # documentation states that 4 cpus is optimal, more yields diminishing returns
        [($encoding eq '64' ? '--phred64-quals' : '--phred33-quals'), '', 0],
        ['-N', $N, 0],
        ['-L', $L, 0],
        [$index_path, '', 0]
    ];

    my ($output_bam) = @$fastq;
    $output_bam =~ s/\.gz$//; # mdb added 7/29/16 -- remove trailing ".fastq" for cutadapt
    $output_bam =~ s/\.fastq$//; # mdb added 7/29/16 -- remove trailing ".fastq" for cutadapt
    $output_bam =~ s/\.fq$//;    # mdb added 8/8/16  -- remove trailing ".fq" for cutadapt

    if ($read_type eq 'paired') {
        $output_bam .= '_bismark_bt2_pe.bam';
        push @$args, ['-1', $fastq->[0], 0];
        push @$args, ['-2', $fastq->[1], 0];
    }
    else { # single-ended
        $output_bam .= '_bismark_bt2.bam';
        push @$args, ['', join(' ', @$fastq), 0];
    }
warn $cmd, $args;
warn $output_bam;
    return {
        cmd         => $cmd,
        args        => $args,
        inputs      => [
            @$fastq,
            @{$self->index},
        ],
        outputs     => [
            $output_bam
        ],
        description => 'Aligning (Bismark) ' . fastq_description($fastq, $read_type)
    };
}

1;