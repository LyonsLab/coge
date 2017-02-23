package CoGe::Builder::Alignment::GSNAP;

use Moose;
extends 'CoGe::Builder::Alignment::Aligner';

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage qw(get_genome_file get_genome_cache_path);
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my %opts = @_;
    my $fasta = $opts{fasta_file}; # reheadered fasta file
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::MissingField->throw(message => 'Missing fastq');
    }

    my $doSeparately = $self->params->{chipseq_params};

    # Generate index
    $self->add(
        $self->gmap_index($fasta)
    );
    $self->index([$self->previous_output->[0]]);

    # Generate sam file
    my @sam_files;
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            $self->add(
                $self->gsnap_alignment([ $file ])
            );
            push @sam_files, $self->previous_output;
        }
    }
    else { # standard GSNAP run (all fastq's at once)
        $self->add(
            $self->gsnap_alignment($fastq)
        );
        push @sam_files, $self->previous_output;
    }

    # Add one or more sam-to-bam tasks
    foreach my $file (@sam_files) {
        # Filter sam file
        $self->add(
            $self->filter_sam($file)
        );

        # Convert sam file to bam
        $self->add(
            $self->sam_to_bam($self->previous_output)
        );
        push @{$self->bam}, $self->previous_output;
    }
}

sub gmap_index {
    my $self = shift;
    my $fasta = shift;

    my $gid = $self->request->genome->id;
    my $name = to_filename($fasta);
    my $GMAP_CACHE_DIR = catdir(get_genome_cache_path($gid), "gmap_index");

    return {
        cmd => get_command_path('GMAP_BUILD'),
        args => [
            ["-D", ".", 0],
            ["-d", $name . "-index", 0],
            ["", $fasta, 1]
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            [catdir($GMAP_CACHE_DIR, $name . "-index"), 1]
        ],
        description => "Indexing genome sequence with GMAP"
    };
}

sub gsnap_alignment {
    my $self = shift;
    my $fastq = shift;
    my $gmap = $self->index->[0];

    my $read_params = $self->params->{read_params} // {};
    my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{read_params} // {};
    my $gapmode = $alignment_params->{'--gap-mode'} // "none";
    my $Q = $alignment_params->{'-Q'} // 1;
    my $n = $alignment_params->{'-n'} // 5;
    my $N = $alignment_params->{'-N'} // 1;
    my $nofails = $alignment_params->{'--nofails'} // 1;
    my $max_mismatches = $alignment_params->{'--max-mismatches'};

    my ($first_fastq) = @$fastq;
    my $output_file = basename($first_fastq) . '.sam';

    my $index_name = basename($gmap);

    my $cmd = get_command_path('GSNAP');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ["-D", ".", 0],
        ["-d", $index_name, 0],
        ["--nthreads=" . $self->NUM_CPUS, '', 0],
        ["-n", $n, 0],
        ["-N", $N, 0],
        ["--format=sam", '', 0],
        ["--gmap-mode=$gapmode", '', 1],
        ["--batch=5", '', 0]
    ];

    push @$args, ["-Q", "", 0] if $Q;
    push @$args, ["--nofails", "", 1] if $nofails;
    push @$args, ["--max-mismatches=$max_mismatches", "", 0] if $max_mismatches;
    push @$args, ['--force-single-end', '', 0] if ($read_type eq 'single');

    # Sort fastq files in case of paired-end reads,
    # see http://research-pub.gene.com/gmap/src/README
    foreach (sort @$fastq) {
        push @$args, ["", $_, 1];
    }

    push @$args, [">", $output_file, 1];

    my $desc = (@$fastq > 2 ? @$fastq . ' files' : join(', ', map { to_filename_base($_) } @$fastq));

    return {
        cmd => $cmd,
# mdb removed 2/2/15 -- fails on zero-length validation input
#        options => {
#            "allow-zero-length" => JSON::false,
#        },
        args => $args,
        inputs => [
            @$fastq,
            [$gmap, 1]
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Aligning $desc with GSNAP"
    };
}

sub filter_sam {
    my $self = shift;
    my $samfile = shift;

    my $filename = basename($samfile);

    my $cmd = catfile($self->conf->{SCRIPTDIR}, "filter_sam.pl");

    return {
        cmd => "$cmd $filename $filename.processed",
        args => [],
        inputs => [
            $samfile
        ],
        outputs => [
            catfile($self->staging_dir, $filename . ".processed")
        ],
        description => "Filtering SAM file"
    };
}

1;