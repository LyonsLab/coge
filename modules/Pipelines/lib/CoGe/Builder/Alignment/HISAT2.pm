package CoGe::Builder::Alignment::HISAT2;

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
    my %opts = @_;
    my $fasta = $opts{fasta_file}; # reheadered fasta file
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $doSeparately = $self->params->{chipseq_params};

    # Add index task
    $self->add(
        $self->hisat2_index($fasta)
    );

    # Add one or more alignment tasks
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            $self->add(
                $self->hisat2_alignment([$file])
            );
            push @sam_files, $self->previous_output;
        }
    }
    else { # standard Bowtie run (all fastq's at once)
        $self->add(
            $self->hisat2_alignment($fastq)
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

sub hisat2_index {
    my $self  = shift;
    my $fasta = shift;

    my $gid = $self->request->genome->id;
    my $cache_dir = catdir(get_genome_cache_path($gid), "hisat2_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');

    my $done_file = "$name.done";

    my $cmd = 'nice ' . get_command_path('HISAT2_BUILD', 'hisat2-build') . " -p 32 $fasta $name && touch $done_file";

    $self->index([
        $name . ".1.ht2",
        $name . ".2.ht2",
        $name . ".3.ht2",
        $name . ".4.ht2",
        $name . ".5.ht2",
        $name . ".6.ht2",
        $name . ".7.ht2",
        $name . ".8.ht2",
    ]);

    return {
        cmd => $cmd,
        args => [],
        inputs => [ $fasta ],
        outputs => [
            @{$self->index},
            $done_file
        ],
        description => "Indexing genome sequence with HISAT2"
    };
}

sub hisat2_alignment {
    my $self  = shift;
    my $fastq = shift;

    my $gid         = $self->request->genome->id;
    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;
    my $read_type   = $read_params->{read_type} // 'single';

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.sam';

	my $args = [
		['-p', $self->NUM_CPUS, 0],
        ['--dta-cufflinks', '', 0], # mdb added 8/9/16 for Cufflinks error "BAM record error: found spliced alignment without XS attribute"
		['-x', catfile(get_genome_cache_path($gid), 'hisat2_index', 'genome.reheader'), 0],
		['-S', $output_file,    0]
    ];

    if ($encoding eq '64') {
        push $args, ['--phred64', '', 0];
    }
    else {
        push $args, ['--phred33', '', 0];
    }

    if ($read_type eq 'single') {
    	push $args, ['-U', join(',', @$fastq), 0];
    }
    else {
		my ($m1, $m2) = detect_paired_end($fastq);
    	push $args, ['-1', join(',', sort @$m1), 0];
    	push $args, ['-2', join(',', sort @$m2), 0];
	}

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

	return {
        cmd => 'nice ' . get_command_path('HISAT2'),
        args => $args,
        inputs => [
            @$fastq,
            @{$self->index}
        ],
        outputs => [ catfile($self->staging_dir, $output_file) ],
        description => "Aligning $desc with HISAT2"
	};
}

1;