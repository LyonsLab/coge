package CoGe::Builder::Alignment::Aligner;

use Moose;
extends 'CoGe::Builder::Buildable';

use Switch;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use Clone qw(clone);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Accessory::Utils qw(is_fastq_file to_filename detect_paired_end to_filename_base to_filename_without_extension);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Core::Metadata qw(to_annotations);
use CoGe::Builder::CommonTasks;
use CoGe::Builder::Trimming::Trimmer;
use CoGe::Builder::Alignment::HISAT2;
use CoGe::Builder::Alignment::GSNAP;
use CoGe::Builder::Alignment::BWA;
use CoGe::Builder::Alignment::Bowtie;
use CoGe::Builder::Alignment::Tophat;
use CoGe::Builder::Alignment::Bismark;
use CoGe::Builder::Alignment::BWAmeth;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

# Outputs
has index   => (is => 'rw', isa => 'ArrayRef', default => sub { [] }); # index files
has raw_bam => (is => 'rw', isa => 'ArrayRef', default => sub { [] }); # unprocessed bam files (straight from the aligner!)
has bam     => (is => 'rw', isa => 'ArrayRef', default => sub { [] }); # processed bam files

sub build {
    my $self = shift;
    my %opts = @_;
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    # Validate inputs not already checked in Request
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }

    my $gid = $self->request->genome->id;
    my $trimming_params  = $self->params->{trimming_params};

# mdb removed 11/6/15 COGE-673
#    # Check multiple files (if more than one file then all should be FASTQ)
#    my $numFastq = 0;
#    foreach (@$input_files) {
#        $numFastq++ if (is_fastq_file($_));
#    }
#    if ($numFastq > 0 and $numFastq != @$input_files) {
#        my $error = 'Unsupported combination of file types';
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }
#    if ($numFastq == 0 and @$input_files > 1) {
#        my $error = 'Too many files';
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }

    #
    # Build workflow
    #

    # Decompress and validate the input files
    my @decompressed;
    foreach my $input_file (@$fastq) {
        if ($input_file =~ /\.gz$/) {
            $self->add_task(
                $self->create_gunzip( $input_file )
            );
            $input_file =~ $self->previous_output(1);

            $self->add_task_chain(
                $self->validate_fastq($input_file)
            );
        }
        else {
            $self->add_task(
                $self->validate_fastq($input_file)
            );
        }

        push @decompressed, $input_file;
    }

    # Trim the fastq input files
    my @trimmed;
    if ($trimming_params) {
        my $trimmer = CoGe::Builder::Trimming::Trimmer->new($self);
        $trimmer->build(data_files => \@decompressed);
        @trimmed = @{$trimmer->fastq};
    }
    else { # no trimming
        @trimmed = @decompressed;
    }

    # Reheader the fasta file
    $self->add_task(
        $self->reheader_fasta($gid)
    );

    # Index the fasta file
    $self->add_task(
        $self->index_fasta($self->previous_output)
    );

    my $aligner;
    switch( $self->_get_aligner() ) {
        case 'hisat2'  { $aligner = CoGe::Builder::Alignment::HISAT2->new($self) }
        case 'bowtie2' { $aligner = CoGe::Builder::Alignment::Bowtie->new($self) }
        case 'tophat'  { $aligner = CoGe::Builder::Alignment::Tophat->new($self) }
        case 'bismark' { $aligner = CoGe::Builder::Alignment::Bismark->new($self) }
        case 'bwameth' { $aligner = CoGe::Builder::Alignment::BWAmeth->new($self) }
        case 'bwa'     { $aligner = CoGe::Builder::Alignment::BWA->new($self) }
        case 'gsnap'   { $aligner = CoGe::Builder::Alignment::GSNAP->new($self) }
        default {
            CoGe::Exception::Generic->throw(message => 'Invalid aligner');
        }
    }
    $aligner->build(data_files => \@trimmed);
    push @{$self->raw_bam}, @{$aligner->bam};

    foreach my $bam_file (@{$self->raw_bam}) {
        # Sort and index the bam output file(s)
        $self->add_task(
            $self->sort_bam($bam_file)
        );
        my $sorted_bam_file = $self->previous_output;
        push @{$self->bam}, $sorted_bam_file;

        $self->add_task(
            $self->index_bam($sorted_bam_file)
        );

        # Get custom metadata to add to experiment #TODO migrate to metadata file
        my $annotations = $self->generate_additional_metadata();

        # Add bam filename to experiment name for ChIP-seq pipeline
        my $md = clone($metadata);
        if (@{$self->raw_bam} > 1) {
            $md->{name} .= ' (' . to_filename_base($sorted_bam_file) . ')';
        }

        # Load alignment
        $self->add_task(
            $self->load_bam(
                metadata => $md,
                annotations => $annotations,
                bam_file => $sorted_bam_file
            )
        );
    }
}

sub validate_fastq {
    my $self = shift;
    my $fastq = shift;

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "validate_fastq.pl"),
        script => undef,
        args => [
            ["", $fastq, 0] # mdb changed 3/1/16 from 1 to 0, COGE-707
        ],
        inputs => [
            $fastq
        ],
        outputs => [
            "$fastq.validated"
        ],
        description => "Validating " . basename($fastq)
    };
}

sub bowtie2_index { # shared between Bowtie and Tophat
    my $self = shift;

    my $gid = $self->request->genome->id;
    my $fasta = get_genome_file($gid);
    my $cache_dir = catdir(get_genome_cache_path($gid), "bowtie_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');

    my $cmd = get_command_path('BOWTIE_BUILD', 'bowtie2-build');

    $self->index([
        $name . ".1.bt2",
        $name . ".2.bt2",
        $name . ".3.bt2",
        $name . ".4.bt2",
        $name . ".rev.1.bt2",
        $name . ".rev.2.bt2"
    ]);

    return catdir($cache_dir, $name), {
        cmd => $cmd,
        args => [
            ["", $fasta, 1],
            ["", $name, 0],
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            @{$self->index}
        ],
        description => "Indexing genome sequence with Bowtie"
    };
}

sub generate_additional_metadata { #TODO redo arg capture in a more automated fashion
    my $self = shift;
    my $read_params      = $self->params->{read_params};
    my $trimming_params  = $self->params->{trimming_params};
    my $alignment_params = $self->params->{alignment_params};

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};

    if ($trimming_params && $trimming_params->{trimmer}) {
        if ($trimming_params->{trimmer} eq 'cutadapt') {
            push @annotations, 'note|cutadapt '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '-m'));
        }
        elsif ($trimming_params->{trimmer} eq 'trimgalore') {
            push @annotations, 'note|trimgalore '. join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--length', '-a'));
        }
    }

    switch( $self->_get_aligner() ) {
        case 'hisat2'  {
            push @annotations, qq{note|hisat2_build};
            my $params = ($read_params->{encoding} eq '64' ? '--phred64' : '--phred33');
            push @annotations, 'note|hisat2 ' . $params;
        }
        case 'bowtie2' {
            my $rg = $alignment_params->{'--rg-id'};
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|bowtie2 ' . $alignment_params->{'presets'} . ($rg ? " --rg-id $rg" : '');
        }
        case 'tophat'  {
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|tophat ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-g'));
        }
        case 'bismark' {
            push @annotations, qq{note|bismark_genome_preparation};
            push @annotations, 'note|bismark ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-L'));
        }
        case 'bwameth' {
            push @annotations, qq{note|bwameth index};
            push @annotations, 'note|bwameth (default options)';
        }
        case 'bwa'     {
            my $M = $alignment_params->{'-M'};
            my $R = $alignment_params->{'-R'};
            my $args_str = ($M ? '-M' : '') . ($R ? " -R $R" : '');
            push @annotations, qq{note|bwa index};
            push @annotations, 'note|bwa mem ' . ($args_str ? $args_str : ' (default options)');
        }
        case 'gsnap'   {
            push @annotations, qq{note|gmap_build};
            push @annotations, 'note|gsnap ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-n', '-Q', '--gap-mode', '--nofails'));
        }
    }

    return \@annotations;
}

sub _get_aligner {
    my $self = shift;
    my $alignment_params = $self->params->{alignment_params};

    if ($alignment_params && $alignment_params->{tool}) {
        return lc($alignment_params->{tool});
    }

    return 'gsnap'; # default aligner if not specified
}

__PACKAGE__->meta->make_immutable;

1;
