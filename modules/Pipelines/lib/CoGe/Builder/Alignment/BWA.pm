package CoGe::Builder::Alignment::BWA;

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
    my $self  = shift;
    my %opts  = @_;
    my $fasta = $opts{fasta_file}; # reheadered fasta file
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $doSeparately = $self->params->{chipseq_params};

    # Add index task
    $self->add(
        $self->bwa_index($fasta)
    );

    # Add one or more alignment tasks
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            $self->add(
                $self->bwa_alignment([$file])
            );
            push @{$self->bam}, $self->previous_output;
        }
    }
    else {
        my $read_params = $self->params->{read_params} // {};
        my $read_type   = $read_params->{read_type} // 'single';
        if ($read_type eq 'paired' && scalar @$fastq > 2) {
            my @bams;
            for (my $i = 0; $i < scalar @$fastq; $i += 2) {
                $self->add(
                    $self->bwa_alignment([$fastq->[$i], $fastq->[$i + 1]])
                );
                $self->add(
                    $self->sort_bam($self->previous_output)
                );
                push @bams, $self->previous_output;
            }
            $self->add(
                $self->merge_bams(\@bams)
            );
            push @{$self->bam}, $self->previous_output;
        }
        else {
            $self->add(
                $self->bwa_alignment($fastq)
            );
            push @{$self->bam}, $self->previous_output;
        }
    }
}

sub bwa_index {
    my $self  = shift;
    my $fasta = shift;

    my $gid = $self->request->genome->id;
    my $cache_dir = catdir(get_genome_cache_path($gid), "bwa_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');
    my $done_file = "$name.done";

    $self->index([
        $name . ".amb",
        $name . ".ann",
        $name . ".bwt",
        $name . ".pac",
        $name . ".sa",
        $done_file
    ]);

    return {
        cmd => 'nice ' . get_command_path('BWA', 'bwa') . ' index',
        args => [
            [ '-a', 'bwtsw', 0 ],
            [ '',   $fasta,  1 ],
            [ '-p', $name,   0 ],
            [ '&&', "touch $done_file", 0 ],
        ],
        inputs => [ $fasta ],
        outputs => [
            @{$self->index},
            $done_file
        ],
        description => "Indexing genome sequence with BWA"
    };
}

sub bwa_alignment {
    my $self  = shift;
    my $fastq = shift;

    my $gid = $self->request->genome->id;

    my $read_params = $self->params->{read_params} // {};
    my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{alignment_params} // {};
    my $M = $alignment_params->{'-M'} // 0;
    my $R = $alignment_params->{'-R'} // '';

    my $CPU = $self->NUM_CPUS;

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.bam';

    my $index_path = catfile(get_genome_cache_path($gid), 'bwa_index', 'genome.reheader');

    my $samtools = get_command_path('SAMTOOLS');

    my @args;
    push @args, ['-M', '', 0] if $M;
    if ($R) {
        $R =~ s/\t/\\t/g; # escape tabs
        push @args, ['-R', shell_quote($R), 0];
    }

    # Build input file list -- decompress bz2 files on-the-fly since bwa only supports decompressed/gzipped files (https://sourceforge.net/p/bio-bwa/mailman/bio-bwa-help/thread/512E3D0C.1030807@bcgsc.ca/)
    my $input_str = join(' ',
        map { (is_bzipped2($_) ? qq['<bunzip2 -c $_'] : $_) }
        sort @$fastq
    );

    push @args, (
        ['-t', $CPU,        0],
        ['',   $index_path, 0],
        ['',   $input_str,  0],
        ["| $samtools view -uSh -\@ $CPU | $samtools sort -\@ $CPU >", $output_file, 1] # convert SAM to BAM and sort on-the-fly for speed
    );

	return {
        cmd => 'nice ' . get_command_path('BWA', 'bwa') . ' mem',
        args => \@args,
        inputs => [
            @$fastq,
            @{$self->index}
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => 'Aligning (BWA-MEM) ' . fastq_description($fastq, $read_type)
	};
}

1;