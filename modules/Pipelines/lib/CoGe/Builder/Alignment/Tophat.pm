package CoGe::Builder::Alignment::Tophat;

use Moose;
extends 'CoGe::Builder::Alignment::Aligner';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my %opts = @_;
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $genome = $self->request->genome;
    my $doSeparately = $self->params->{chipseq_params};

    # Generate gff if genome annotated
    my $gff_file;
    if ( $genome->has_gene_features ) {
        $self->add(
            $self->create_gff( #FIXME simplify this
                gid => $genome->id,
                output_file => get_gff_cache_path(
                    gid => $genome->id,
                    genome_name => sanitize_name($genome->organism->name),
                    output_type => 'gff',
                    params => {}
                )
            )
        );
        $gff_file = $self->previous_output;
    }

    # Add index task
    $self->add(
        $self->bowtie2_index()
    );

    # Add one or more alignment tasks
    if ($doSeparately) { # ChIP-seq pipeline (align each fastq individually)
        foreach my $file (@$fastq) {
            $self->add(
                $self->tophat_alignment([$file], $gff_file)
            );
            push @{$self->bam_files}, $self->previous_output;
        }
    }
    else { # standard Tophat run (all fastq's at once)
        $self->add(
            $self->tophat_alignment($fastq, $gff_file)
        );
        push @{$self->bam_files}, $self->previous_output;
    }
}

sub create_tophat_job {
    my $self  = shift;
    my $fastq = shift;
    my $gff   = shift; # optional

    my $gid         = $self->request->genome->id;
    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;
    my $read_type   = $read_params->{read_type} // 'single';

    my $alignment_params = $self->params->{alignment_params} // {};
    my $g = $alignment_params->{'-g'} // 1;

    my $index_name  = catfile(get_genome_cache_path($gid), 'bowtie_index', 'genome.reheader');

    # Setup input dependencies
    my $inputs = [
        @$fastq,
        @{$self->index}
    ];
    push @$inputs, $gff if $gff;

    my ($first_fastq) = @$fastq;
    my $output_file = to_filename_without_extension($first_fastq) . '.bam';

    # Build up command/arguments string
    my $cmd = get_command_path('TOPHAT');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $arg_str;
    $arg_str .= $cmd . ' ';
    $arg_str .= "-G $gff " if ($gff);
    $arg_str .= "--phred64_quals " if ($encoding eq '64');
    $arg_str .= "-o . -g $g -p 32 $index_name ";

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'tophat.pl'), # this script was created because JEX can't handle TopHat's paired-end argument syntax
        script => undef,
        args => [
            ['-read_type', $read_type, 0],
            ['-cmd_args', shell_quote($arg_str), 0],
            ['-output', $output_file, 0],
            ['-files', shell_quote(join(',', @$fastq)), 0]
        ],
        inputs => $inputs,
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => "Aligning sequences with TopHat"
    };
}

1;