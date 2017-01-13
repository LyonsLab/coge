package CoGe::Builder::Trimming::Trimmer;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);

use CoGe::Accessory::Utils qw(detect_paired_end);
use CoGe::Builder::Trimming::Cutadapt;
use CoGe::Builder::Trimming::TrimGalore;
use CoGe::Exception::Generic;

# Outputs
has fastq => (is => 'rw', isa => 'ArrayRef', default => sub { [] }); # trimmed fastq files

sub build {
    my $self  = shift;
    my %opts = @_;
    my $fastq = $opts{data_files}; # array ref of fastq files
    unless ($fastq && @$fastq) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $read_params      = $self->params->{read_params};
    my $trimming_params  = $self->params->{trimming_params};

    my ($fastq1, $fastq2); # array refs
    if ($read_params->{read_type} eq 'paired') {
        # Separate files based on last occurrence of _R1 or _R2 in filename
        my ($m1, $m2) = detect_paired_end($fastq);
        unless (@$m1 and @$m2 and @$m1 == @$m2) {
            CoGe::Exception::Generic->throw(message => 'Mispaired FASTQ files', details => Dumper { m1 => $m1, m2 => $m2 });
        }
        $fastq1 = $m1;
        $fastq2 = $m2;
    }
    else { # default to single-ended
        $fastq1 = $fastq;
    }

    my $trimmer;
    if ($trimming_params->{trimmer} eq 'cutadapt') {
        $trimmer = CoGe::Builder::Trimming::Cutadapt->new($self);
    }
    elsif ($trimming_params->{trimmer} eq 'trimgalore') {
        $trimmer = CoGe::Builder::Trimming::TrimGalore->new($self);
    }
    else {
        CoGe::Exception::Generic->throw(message => 'Unrecognized trimmer');
    }

    $trimmer->build(data_files => [$fastq1, $fastq2]);
    push @{$self->fastq}, @{$trimmer->fastq};
}

1;