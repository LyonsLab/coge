package CoGe::Builder::Protein::ChIPseq;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

my $HOMER_DIR;
sub BUILD { # called immediately after constructor
    my $self = shift;
    $HOMER_DIR = $self->conf->{HOMER_DIR};
    unless ($HOMER_DIR) {
        CoGe::Exception::Generic->throw(message => 'Missing HOMER_DIR in config file');
    }
}

sub build {
    my $self = shift;
    my $bam_files = shift;
    unless ($bam_files && @$bam_files == 3) { # use input experiment's bam file (for MeasureExpression)
        CoGe::Exception::Generic->throw(message => 'Missing bam inputs');
    }

    unless ($self->params->{metadata}) { # use input experiment's metadata (for MeasureExpression)
        my $experiment = self->request->experiment;
        $self->params->{metadata} = { # could almost use experiment->to_hash here except for source_name
            name       => $experiment->name,
            version    => $experiment->version,
            source     => $experiment->source->name,
            restricted => $experiment->restricted
        };
    }

    my $genome = $self->request->genome;

    # Set metadata for the pipeline being used
    my $annotations = generate_additional_metadata();
    my @annotations2 = CoGe::Core::Metadata::to_annotations($additional_metadata);
    push @$annotations, @annotations2;

    my $chipseq_params = $self->params->{chipseq_params};
    unless (chipseq_params->{input}) {
        CoGe::Exception::Generic->throw(message => 'Missing input designation');
    }

    # Determine which bam file corresponds to the input vs. the replicates
    my $input_base = to_filename_base($chipseq_params->{input});
    my ($input_file, @replicates);
    foreach my $file (@$bam_files) {
        my ($basename) = to_filename_base($file);
        if (index($basename, $input_base) != -1) {
            $input_file = $file;
        }
        else {
            push @replicates, $file;
        }
    }
    unless ($input_file) {
        CoGe::Exception::Generic->throw( message => "Unable to detect input, base=$input_base", details => Dumper $input_files);
    }

    #
    # Build the workflow
    #
    foreach my $bam_file (@$bam_files) {
        $self->add_task(
            $self->bamToBed($bam_file)
        );

        $self->add_task(
            $self->homer_makeTagDirectory(
                bed_file => $self->previous_output,
                gid => $genome->id
            )
        );
    }
    
    foreach my $replicate (@replicates) {
        my ($input_tag) = to_filename_base($input_file);
        my ($replicate_tag) = to_filename_base($replicate);

        $self->add_task(
            $self->homer_findPeaks(
                input_dir => catdir($self->staging_dir, $input_tag),
                replicate_dir => catdir($self->staging_dir, $replicate_tag)
            )
        );

        $self->add_task(
            $self->convert_homer_to_csv($self->previous_output)
        );

        my $md = clone($metadata);
        $md->{name} .= " ($input_tag vs. $replicate_tag) (ChIP-seq)";
        push @{$md->{tags}}, 'ChIP-seq';

        $self->add_task(
            $self->load_experiment(
                metadata    => $md,
                gid         => $genome->id,
                input_file  => $self->previous_outputs,
                name        => $replicate_tag,
                normalize   => 'percentage',
                annotations => $annotations
            )
        );
    }
}

sub bamToBed {
    my $self = shift;
    my $bam_file = shift;

    my $cmd = $self->conf->{BAMTOBED} || 'bamToBed';
    
    my $name = to_filename_base($bam_file);
    my $bed_file = catfile($self->staging_dir, $name . '.bed');
    my $done_file = $bed_file . '.done';
    
    return {
        cmd => "$cmd -i $bam_file > $bed_file ; touch $done_file",
        args => [],
        inputs => [
            $bam_file
        ],
        outputs => [
            $bed_file,
            $done_file
        ],
        description => "Converting $name BAM file to BED format"
    };
}

sub homer_makeTagDirectory {
    my $self = shift;
    my %opts = @_;
    my $bed_file = $opts{bed_file};
    my $gid      = $opts{gid};

    my $params = $self->params->{chipseq_params};
    my $size   = $params->{'-size'} // 250;

    my $tag_name = to_filename_base($bed_file);
    
    my $fasta = catfile(get_genome_cache_path($gid), 'genome.faa.reheader.faa'); #TODO move into function in Storage.pm
    
    return {
        cmd => catfile($HOMER_DIR, 'makeTagDirectory'),
        script => undef,
        args => [
            ['', $tag_name, 0],
            ['', $bed_file, 0],
            ['-fragLength', $size, 0],
            ['-format', 'bed', 0],
            ['-genome', $fasta, 0],
            ['-checkGC', '', 0]
        ],
        inputs => [
            $bed_file,
            qq[$bed_file.done],
            $fasta
        ],
        outputs => [
            [catfile($self->staging_dir, $tag_name), 1],
            catfile($self->staging_dir, $tag_name, 'tagInfo.txt')
        ],
        description => "Creating tag directory '$tag_name' using Homer"
    };
}

sub homer_findPeaks {
    my $self = shift;
    my %opts = @_;
    my $replicate_dir = $opts{replicate_dir};
    my $input_dir     = $opts{input_dir};
    my $staging_dir   = $opts{staging_dir};

    my $params = $self->params->{chipseq_params};
    my $size   = $params->{'-size'}  // 250;
    my $gsize  = $params->{'-gsize'} // 3000000000;
    my $norm   = $params->{'-norm'}  // 1e8;
    my $fdr    = $params->{'-fdr'}   // 0.01;
    my $F      = $params->{'-F'}     // 3;

    my ($replicate_tag) = to_filename_base($replicate_dir);
    
    my $output_file = "homer_peaks_$replicate_tag.txt";
    
    return {
        cmd => catfile($HOMER_DIR}, 'findPeaks'),
        script => undef,
        args => [
            ['',       $replicate_dir, 0],
            ['-i',     $input_dir,     0],
            ['-style', 'factor',       0],
            ['-o',     $output_file,   0],
            ['-size',  $size,          0],
            ['-gsize', $gsize,         0],
            ['-norm',  $norm,          0],
            ['-fdr',   $fdr,           0],
            ['-F',     $F,             0]
        ],
        inputs => [
            [$replicate_dir, 1],
            [$input_dir, 1]
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Performing ChIP-seq analysis on $replicate_tag using Homer"
    };
}

sub convert_homer_to_csv {
    my $self = shift;
    my $input_file = shift;

    my $cmd = catfile($self->conf->{SCRIPTDIR}, 'chipseq', 'homer_peaks_to_csv.pl');
    
    my $name = to_filename_without_extension($input_file);
    my $output_file = catfile($self->staging_dir, $name . '.csv');
    my $done_file = $output_file . '.done';
    
    return {
        cmd => "$cmd $input_file > $output_file ; touch $done_file",
        script => undef,
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            $done_file
        ],
        description => "Converting $name to CSV format"
    };
}

sub generate_additional_metadata {
    my $self = shift;
    my $chipseq_params = $self->params->{chipseq_params};
    $chipseq_params->{'-fragLength'} = $chipseq_params->{'-size'}; # kludge b/c "size" is used for "fragLength" argument

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExperiment||note|Generated by CoGe's NGS Analysis Pipeline};
    push @annotations, 'note|makeTagDirectory ' . join(' ', map { $_.' '.$chipseq_params->{$_} } ('-fragLength', '-checkGC'));
    push @annotations, 'note|findPeaks ' . join(' ', map { $_.' '.$chipseq_params->{$_} } ('-size', '-gsize', '-norm', '-fdr', '-F'));

    return \@annotations;
}

1;