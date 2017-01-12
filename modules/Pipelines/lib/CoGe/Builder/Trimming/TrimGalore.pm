package CoGe::Builder::Trimming::TrimGalore;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use Data::Dumper qw(Dumper);
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my $fastq1 = shift;
    my $fastq2 = shift; # undef for single-ended
    unless ($fastq1 && @$fastq1) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    if ($read_type eq 'single') { # single-ended
        # Create cutadapt task for each file
        foreach my $file (@$fastq1) {
            $self->add_task(
                $self->trimgalore($file)
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create cutadapt task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my $file1 = shift @$fastq1;
            my $file2 = shift @$fastq2;
            $self->add_task(
                $self->trimgalore([ $file1, $file2 ])
            );
            push @{$self->fastq}, $self->previous_outputs;
        }
    }
}

sub create_trimgalore_job {
    my $self = shift;
    my $fastq = shift; # for single fastq file (backwards compatibility) or two paired-end fastq files (new functionality)

    my $trimming_params = $self->params->{trimming_params} // {};
    my $q = $trimming_params->{'-q'} // 20;
    my $length = $trimming_params->{'-length'} // 20;
    my $a = $trimming_params->{'-a'};
    my $read_params = $self->params->{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');

    my $name = join(', ', map { basename($_) } @$fastq);

    # Build up command/arguments string
    my $cmd = $self->conf->{TRIMGALORE} || 'trim_galore';
    $cmd = 'nice ' . $cmd; # run at lower priority

    # Create staging dir
    $cmd = "mkdir -p $self->staging_dir && " . $cmd;

    my $args = [
        ['--output_dir', $self->staging_dir, 0],
        ['-q', $q, 0],
        ['--length', $length, 0],
    ];

    my $phred = ($encoding == 64 ? '--phred64' : '--phred33');
    push @$args, [$phred, '', 0];

    if (defined $a) {
        push @$args, ['-a', '"'.$a.'"', 0];
    }

    my @outputs;
    if ($read_type eq 'paired') {
        push @$args, ['--paired', join(' ', @$fastq), 0];

        my ($r1, $r2) = @$fastq;
        @outputs = ( catfile($self->staging_dir, to_filename_without_extension($r1) . '_val_1.fq'),
                     catfile($self->staging_dir, to_filename_without_extension($r2) . '_val_2.fq') );
    }
    else { # single
        my ($file) = @$fastq;
        push @$args, ['', $file, 0];
        @outputs = ( catfile($self->staging_dir, to_filename_without_extension($file) . '_trimmed.fq') );
    }

    # kludge to fix JEX sequencing
    foreach (@$args) {
        $cmd .= ' ' . $_->[0] . ' ' . $_->[1];
    }
    my @done_files = map { $_ . '.done' } @outputs;
    foreach (@done_files) {
        $cmd .= " && touch $_";
    }

    return {
        cmd => $cmd,
        args => [],
        inputs => [
            @$fastq
        ],
        outputs => \@outputs,
        description => "Trimming (trimgalore) $name"
    };
}

1;