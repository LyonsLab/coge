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
    my %opts = @_;
    my $fasta = $opts{fasta_file}; # reheadered fasta file
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    # Add index task
    $self->add(
        $self->bwameth_index($fasta)
    );

    $self->add_to_previous(
        $self->bwameth_alignment($fastq)
    );
    push @{$self->bam}, $self->previous_output;
}

sub bwameth_index {
    my $self = shift;
    my $fasta = shift;

    my $name = to_filename($fasta);

    my $gid = $self->request->genome->id;
    my $cache_dir = catdir(get_genome_cache_path($gid), "bwameth_index");
	make_path($cache_dir) unless (-d $cache_dir);

    my $done_file = 'bwameth_index.done';

    my $cmd = ($self->conf->{BWAMETH} ? $self->conf->{BWAMETH} : 'bwameth') . ' index';

    $cmd = "cd $cache_dir && " .
           "cp $fasta . && " .
           "$cmd $name && " .
           "touch $done_file";

    $self->index([
        catfile($cache_dir, "$name.bwameth.c2t"),
        catfile($cache_dir, "$name.bwameth.c2t.amb"),
        catfile($cache_dir, "$name.bwameth.c2t.ann"),
        catfile($cache_dir, "$name.bwameth.c2t.bwt"),
        catfile($cache_dir, "$name.bwameth.c2t.pac"),
        catfile($cache_dir, "$name.bwameth.c2t.sa")
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
    my $index_path  = catdir(get_genome_cache_path($gid), "bwameth_index");

    my $cmd = $self->conf->{BWAMETH} || 'bwameth';
    $cmd = 'nice ' . $cmd; # run at lower priority

    return {
        cmd         => $cmd,
        args        => [
            ['--reference', catfile($index_path, 'genome.faa.reheader.faa'), 0],
            ['', join(' ', @$fastq), 0],
            ['-t', 8, 0],
            ['-p', 'alignment', 0]
        ],
        inputs      => [
            @$fastq,
            @{$self->index}
        ],
        outputs     => [
            catfile($self->staging_dir, 'alignment.bam')
        ],
        description => "Aligning sequences with bwameth"
    };
}

1;