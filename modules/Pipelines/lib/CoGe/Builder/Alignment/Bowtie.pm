package CoGe::Builder::Alignment::Bowtie;

use Moose;
extends 'CoGe::Builder::Alignment::Aligner';

use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

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

    my $doSeparately = $self->params->{chipseq_params};

    # Add index task
    $self->add(
        $self->bowtie2_index($fasta)
    );

    # Add one or more alignment tasks
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            $self->add(
                $self->bowtie2_alignment([$file])
            );
            push @sam_files, $self->previous_output;
        }
    }
    else { # standard Bowtie run (all fastq's at once)
        $self->add(
            $self->bowtie2_alignment($fastq)
        );
        push @sam_files, $self->previous_output;
    }

    # Add one or more sam-to-bam tasks
    foreach my $file (@sam_files) {
        $self->add(
            $self->sam_to_bam($file)
        );
        push @{$self->bam}, $self->previous_output;
    }
}

sub bowtie2_alignment {
    my $self  = shift;
    my $fastq = shift;

    my $gid         = $self->request->genome->id;
    my $index_name  = catfile(get_genome_cache_path($gid), 'bowtie_index', 'genome.reheader');

    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;
    my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{alignment_params} // {};
    my $presets     = $alignment_params->{presets} // '--sensitive';
    my $read_group  = $alignment_params->{read_group} // '';

    # Setup input dependencies
    my $inputs = [
        @$fastq,
        @{$self->index}
    ];

    # Build up command/arguments string
    my $cmd = $self->conf->{BOWTIE2} || 'bowtie2';
    $cmd = 'nice ' . $cmd; # run at lower priority
    $cmd .= ' -p ' . $self->NUM_CPUS;
    $cmd .= ' ' . shell_quote($presets);
    $cmd .= ' --rg-id ' . shell_quote($read_group) if $read_group;
    $cmd .= ' --phred64' if ($encoding eq '64'); # default is --phred33
    $cmd .= " -x $index_name ";

    if ($read_type eq 'paired') {
        my ($m1, $m2) = detect_paired_end(\@$fastq);
        unless (@$m1 and @$m2) {
            CoGe::Exception::Generic->throw( message => 'Mispaired FASTQ files', details => Dumper { m1 => $m1, m2 => $m2 } );
        }
        $cmd .= '-1 ' . join(',', sort @$m1) . ' -2 ' . join(',', sort @$m2);
    }
    else { # single-ended
        $cmd .= '-U ' . join(',', sort @$fastq);
    }

    my $samtools = get_command_path('SAMTOOLS');
    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.bam';
    $cmd .= " | $samtools view -bS > $output_file"; # convert SAM to BAM on the fly for speed

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    return {
        cmd => $cmd,
        args => [],
        inputs => $inputs,
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Aligning $desc with Bowtie2"
    };
}

1;