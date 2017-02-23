package CoGe::Builder::Alignment::BWA;

use Moose;
extends 'CoGe::Builder::Alignment::Aligner';

use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils qw(to_filename to_filename_base to_filename_without_extension);
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
    else { # standard Bowtie run (all fastq's at once)
        $self->add(
            $self->bwa_alignment($fastq)
        );
        push @{$self->bam}, $self->previous_output;
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
    
    #my $read_params = $self->params->{read_params} // {};
    #my $encoding    = $read_params->{encoding} // 33;
    #my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{alignment_params} // {};
    my $M = $alignment_params->{'-M'} // 0;
    my $R = $alignment_params->{'-R'} // '';

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.bam';

    my $index_path = catfile(get_genome_cache_path($gid), 'bwa_index', 'genome.reheader');

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    my $samtools = get_command_path('SAMTOOLS');

    my @args;
    push @args, ['-M', '', 0] if $M;
    push @args, ['-R', shell_quote($R), 0] if $R;
    push @args, (
        ['-t', $self->NUM_CPUS,         0],
        ['',   $index_path,             0],
        ['',   join(' ', sort @$fastq), 0],
        ["| $samtools view -bS >", $output_file, 1] # convert SAM to BAM on the fly for speed
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