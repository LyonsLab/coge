package CoGe::Builder::Trimming::Trimmomatic;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use Data::Dumper;
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my %opts = @_;
    my ($fastq1, $fastq2) =  @{$opts{data_files}}; # fastq2 is undef for single-ended
    unless ($fastq1 && @$fastq1) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    if ($read_type eq 'single') { # single-ended
        # Create trimmomatic task for each file
        foreach (@$fastq1) {
            $self->add(
                $self->trimmomatic($_),
                qq[$_.done] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create trimmomatic task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my ($f1, $f2) = ($fastq1->[$i], $fastq2->[$i]);
            $self->add(
                $self->trimmomatic([ $fastq1->[$i], $fastq2->[$i] ]),
                [ qq{$f1.done}, qq{$f2.done} ] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, grep { /\.fastq$/ } @{$self->previous_outputs};
        }
    }
}

sub trimmomatic {
    my $self = shift;
    my $fastq = shift; # single fastq file or two paired-end fastq files
    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');

    my $trimming_params = $self->params->{trimming_params} // {};
    my $read_params = $self->params->{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    my @outputs = map { catfile($self->staging_dir, basename(remove_fastq_ext($_) . '.trimmed.fastq' . to_compressed_ext($_))) } @$fastq; #map { catfile($self->staging_dir, to_filename($_) . '.trimmed.fastq') } @$fastq;

    my $cmd = get_command_path('TRIMMOMATIC'); # path to jar file
    $cmd = 'nice java -jar ' . $cmd; # run at lower priority

    my $args;

    if ($read_type eq 'paired') { # paired-end
        # Example:  java -jar trimmomatic-0.36.jar PE SRR1564620_1.fastq SRR1564620_2.fastq -baseout trimmed.fastq MINLEN:36
        push @$args, [ 'PE', '', 0 ];
        push @$args, map { [ '', $_, 1 ] } @$fastq;
        push @$args, [ '-baseout', to_filename($fastq->[0]), 0 ];
    }
    else { # single-ended
        # Example:  java -jar trimmomatic-0.36.jar SE SRR1564620.fastq trimmed MINLEN:36
        push @$args, [ 'SE', $fastq->[0], 1 ];
        push @$args, [ '',   $outputs[0], 1 ];
    }

    foreach ( 'ILLUMINACLIP', 'SLIDINGWINDOW', 'MAXINFO', 'LEADING', 'TRAILING', 'CROP', 'HEADCROP', 'MINLEN' ) {
        my $value = $trimming_params->{$_};
        push @$args, ["$_:$value", '', 0] if ($value);
    }

    push @$args, ( # these must go at end of command-line
        [ '-threads', 8, 0 ],
        [ '-trimlog', 'trimmomatic.log', 0 ],
        [($encoding == 64 ? '-phred64' : '-phred33'), '', 0]
    );

    return {
        cmd => $cmd,
        args => $args,
        inputs => [ @$fastq ],
        outputs => \@outputs,
        description => 'Trimming (trimmomatic) ' . join(', ', map { basename($_) } @$fastq)
    };
}

1;