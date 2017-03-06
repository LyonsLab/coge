package CoGe::Builder::Trimming::Cutadapt;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use Data::Dumper;
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);
use String::ShellQuote qw(shell_quote);

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
        # Create cutadapt task for each file
        foreach (@$fastq1) {
            $self->add(
                $self->cutadapt($_),
                qq[$_.done] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create cutadapt task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my ($f1, $f2) = ($fastq1->[$i], $fastq2->[$i]);
            $self->add(
                $self->cutadapt([ $f1, $f2 ]),
                [ qq{$f1.done}, qq{$f2.done} ] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, @{$self->previous_outputs};
        }
    }
}

sub cutadapt {
    my $self = shift;
    my $fastq = shift; # single fastq file or two paired-end fastq files
    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');

    my $trimming_params = $self->params->{trimming_params} // {};
    my $q = $trimming_params->{'-q'} // 25;
    my $m = $trimming_params->{'-m'} // 17;
    my $read_params = $self->params->{read_params} // {};
    my $encoding  = $read_params->{encoding}  // 33;
    my $read_type = $read_params->{read_type} // 'single';

    my @outputs = map { catfile($self->staging_dir, basename(remove_fastq_ext($_) . '.trimmed.fastq' . to_compressed_ext($_))) } @$fastq;

    my $cmd = get_command_path('CUTADAPT');
    $cmd = 'nice ' . $cmd; # run at lower priority

    my $arg_str  = "$cmd -q $q --quality-base=$encoding -m $m -o $outputs[0]";
       $arg_str .= " -p $outputs[1] " if ($read_type eq 'paired'); # paired-end

    my $args = [
        [ $read_type,            '', 0 ],
        [ $self->staging_dir,    '', 0 ],
        [ shell_quote($arg_str), '', 0 ]
    ];
    push @$args, map { [ '', $_, 1 ] } @$fastq; # relative path

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'cutadapt.pl'), # this script was created because JEX can't handle Cutadapt's paired-end argument syntax
        args => $args,
        inputs => [ @$fastq ],
        outputs => \@outputs,
        description => 'Trimming (cutadapt) ' . join(', ', map { basename($_) } @$fastq)
    };
}

1;