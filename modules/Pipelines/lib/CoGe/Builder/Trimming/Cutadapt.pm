package CoGe::Builder::Trimming::Cutadapt;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use Data::Dumper;
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils qw(to_filename);
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
        # Create cutadapt task for each file
        foreach my $file (@$fastq1) {
            $self->add_task(
                $self->cutadapt($file)
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create cutadapt task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            $self->add_task(
                $self->cutadapt([ $fastq1->[$i], $fastq2->[$i] ])
            );
            push @{$self->fastq}, $self->previous_outputs;
        }
    }
}

sub cutadapt {
    my $self = shift;
    my $fastq = shift; # for single fastq file (backwards compatibility) or two paired-end fastq files (new functionality)

    my $trimming_params = $self->params->{trimming_params} // {};
    my $q = $trimming_params->{'-q'} // 25;
    my $m = $trimming_params->{'-m'} // 17;
    my $read_params = $self->params->{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');

    my $name = join(', ', map { basename($_) } @$fastq);
    my @outputs = map { catfile($self->staging_dir, to_filename($_) . '.trimmed.fastq') } @$fastq;

    # Build up command/arguments string
    my $cmd = get_command_path('CUTADAPT');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $arg_str;
    $arg_str .= $cmd . ' ';
    $arg_str .= "-q $q --quality-base=$encoding -m $m -o $outputs[0] ";
    $arg_str .= "-p $outputs[1] " if ($read_type eq 'paired'); # paired-end

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'cutadapt.pl'), # this script was created because JEX can't handle Cutadapt's paired-end argument syntax
        args => [
            [ $read_type, '', 0 ],
            [ $self->staging_dir, '', 0 ],
            [ '"'.$arg_str.'"', '', 0 ],
            [ '', join(' ', @$fastq), 0 ]
        ],
        inputs => [
            @$fastq
        ],
        outputs => \@outputs,
        description => "Trimming (cutadapt) $name"
    };
}

1;