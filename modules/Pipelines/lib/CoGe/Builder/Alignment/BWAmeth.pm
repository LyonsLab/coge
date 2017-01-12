package CoGe::Builder::Alignment::BWAmeth;

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

    # Add index task
    $self->add_task(
        $self->bwameth_index()
    );

    $self->add_task(
        $self->bwameth_alignment($fastq)
    );
    push @{$self->bam}, $self->previous_output;
}

sub bwameth_index {
    my $self = shift;

    my $gid = $self->request->genome->id;
    my $fasta = get_genome_file($gid);
    my $cache_dir = catdir(get_genome_cache_path($gid), "bwameth_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');

    my $done_file = 'bwameth_index.done';

    my $cmd = ($self->conf->{BWAMETH} ? $self->conf->{BWAMETH} : 'bwameth') . ' index';

    $cmd = "cd $cache_dir && " .
           "cp $fasta . && " .
           "$cmd $name && " .
           "touch $done_file";

    $self->index([
        "$name.bwameth.c2t",
        "$name.bwameth.c2t.amb",
        "$name.bwameth.c2t.ann",
        "$name.bwameth.c2t.bwt",
        "$name.bwameth.c2t.pac",
        "$name.bwameth.c2t.sa"
    ]);

    return {
        cmd => $cmd,
        args => [],
        inputs => [
            $fasta
        ],
        outputs => [
            @{$self->index},
            catfile($cache_dir, $done_file)
        ],
        description => "Indexing genome sequence with bwameth"
    };
}

sub bwameth_alignment {
    my $self  = shift;
    my $fastq = shift;

    my $gid         = $self->request->genome->id;
    my $read_params = $self->params->{read_params} // {};

    my $index_path  = catdir(get_genome_cache_path($gid), "bwameth_index");

    # Setup input dependencies
    my $inputs = [
        @$fastq,
        @{$self->index}
    ];

    # Build command and arguments
    my $cmd = $self->conf->{BWAMETH} || 'bwameth';
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ['--reference', catfile($index_path, 'genome.faa.reheader.faa'), 0],
        ['', join(' ', @$fastq), 0],
        ['-t', 8, 0],
        ['-p', 'alignment', 0]
    ];

    return (
        cmd => $cmd,
        args => $args,
        inputs => $inputs,
        outputs => [
            catfile($self->staging_dir, 'alignment.bam')
        ],
        description => "Aligning sequences with bwameth"
    );
}

1;