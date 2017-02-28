package CoGe::Builder::Trimming::TrimGalore;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use Data::Dumper qw(Dumper);
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my %opts = @_;
    my ($fastq1, $fastq2) = @{$opts{data_files}}; # fastq2 is undef for single-ended
    unless ($fastq1 && @$fastq1) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    if ($read_type eq 'single') { # single-ended
        # Create cutadapt task for each file
        foreach (@$fastq1) {
            $self->add(
                $self->trimgalore($_)
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create cutadapt task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            $self->add(
                $self->trimgalore([ $fastq1->[$i], $fastq2->[$i] ])
            );
            push @{$self->fastq}, $self->previous_outputs;
        }
    }
}

sub trimgalore {
    my $self = shift;
    my $fastq = shift; # single fastq file or two paired-end fastq files
    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');

    my $trimming_params = $self->params->{trimming_params} // {};
    my $q      = $trimming_params->{'-q'} // 20;
    my $length = $trimming_params->{'-length'} // 20;
    my $a      = $trimming_params->{'-a'};
    my $read_params = $self->params->{read_params} // {};
    my $encoding    = $read_params->{encoding} // 33;
    my $read_type   = $read_params->{read_type} // 'single';

    # Build up command/arguments string
    my $cmd = $self->conf->{TRIMGALORE} || 'trim_galore';
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $args = [
        ['--output_dir', $self->staging_dir, 0],
        ['-q',           $q,                 0],
        ['--length',     $length,            0],
        [($encoding == 64 ? '--phred64' : '--phred33'), '', 0]
    ];

    push @$args, ['-a', '"'.$a.'"', 0] if (defined $a);

    my @outputs;
    if ($read_type eq 'paired') {
        push @$args, ['--paired', '', 0];
        push @$args, map { ['', $_, 1] } @$fastq; # relative path

        my ($r1, $r2) = @$fastq;
        @outputs = ( catfile($self->staging_dir, basename(remove_fastq_ext($r1) . '_val_1.fq' . to_compressed_ext($r1))),   #catfile($self->staging_dir, to_filename_without_extension($r1) . '_val_1.fq'),
                     catfile($self->staging_dir, basename(remove_fastq_ext($r2) . '_val_2.fq' . to_compressed_ext($r2))) ); #catfile($self->staging_dir, to_filename_without_extension($r2) . '_val_2.fq') );
    }
    else { # single
        push @$args, map { ['', $_, 1] } @$fastq; # relative path
        @outputs = map { catfile($self->staging_dir, basename(remove_fastq_ext($_) . '_trimmed.fq' . to_compressed_ext($_))) } @$fastq; #map { catfile($self->staging_dir, to_filename_without_extension($_) . '_trimmed.fq') } @$fastq;
    }

#    # kludge to fix JEX sequencing
#    foreach (@$args) {
#        $cmd .= ' ' . $_->[0] . ' ' . $_->[1];
#    }
#    my @done_files = map { $_ . '.done' } @outputs;
#    foreach (@done_files) {
#        $cmd .= " && touch $_";
#    }

    return {
        cmd => $cmd,
        args => $args,
        inputs => [ @$fastq ],
        outputs => \@outputs,
        description => 'Trimming (trimgalore) ' . join(', ', map { basename($_) } @$fastq)
    };
}

1;