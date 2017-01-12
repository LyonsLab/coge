package CoGe::Builder::Alignment::BWA;

use Moose;
extends 'CoGe::Builder::Alignment::Aligner';

use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);

use CoGe::Accessory::Utils qw(to_filename to_filename_base to_filename_without_extension);
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage qw(get_genome_file get_genome_cache_path);
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my $fastq = shift; # array ref of fastq files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $doSeparately = $self->params->{chipseq_params};

    # Add index task
    $self->add_task(
        $self->bwa_index()
    );

    # Add one or more alignment tasks
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            $self->add_task(
                $self->bwa_alignment([$file])
            );
            push @sam_files, $self->previous_output;
        }
    }
    else { # standard Bowtie run (all fastq's at once)
        $self->add_task(
            $self->bwa_alignment($fastq)
        );
        push @sam_files, $self->previous_output;
    }

    # Add one or more sam-to-bam tasks
    foreach my $file (@sam_files) {
        $self->add_task(
            $self->sam_to_bam($file)
        );
        push @{$self->bam}, $self->previous_output;
    }
}

sub bwa_index {
    my $self  = shift;
    my $gid = $self->request->genome->id;
    my $fasta = get_genome_file($gid);

    my $cache_dir = catdir(get_genome_cache_path($gid), "bwa_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');
    my $done_file = "$name.done";

    my $cmd = 'nice ' . get_command_path('BWA', 'bwa') . " index $fasta -p $name && touch $done_file";

    $self->index([
        $name . ".amb",
        $name . ".ann",
        $name . ".bwt",
        $name . ".pac",
        $name . ".sa",
    ]);

    return {
        cmd => $cmd,
        args => [],
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

    my $gid         = $self->request->genome->id;
    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;
    my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{alignment_params} // {};
    my $M = $alignment_params->{'-M'} // 0;
    my $R = $alignment_params->{'-R'} // '';

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.sam';

    my $index_path = catfile(get_genome_cache_path($gid), 'bwa_index', 'genome.reheader');

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    my @args;
    push @args, ['-M', '', 0] if $M;
    push @args, ['-R', shell_quote($R), 0] if $R;
    push @args, (
        ['-t', '32',                    0],
        ['',   $index_path,             0],
        ['',   join(' ', sort @$fastq), 0],
        ['>',  $output_file,            1]
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
        description => "Aligning $desc using BWA-MEM"
	};
}

1;