package CoGe::Builder::Trimming::Trimmer;

use Moose;
extends 'CoGe::Builder::Buildable';

use Switch;
use Data::Dumper qw(Dumper);

use CoGe::Accessory::Utils qw(detect_paired_end);
use CoGe::Builder::Trimming::Cutadapt;
use CoGe::Builder::Trimming::TrimGalore;
use CoGe::Builder::Trimming::Trimmomatic;
use CoGe::Builder::Trimming::BBDuk;
use CoGe::Exception::Generic;

# Settings
has NUM_CPUS => (is => 'ro', isa => 'Int', default => 8); # number of CPUs to use for trimming tasks (when applicable)

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
            CoGe::Exception::Generic->throw(message => 'Mispaired FASTQ files. Filenames must include _R1 and _R2.');
        }
        $fastq1 = $m1;
        $fastq2 = $m2;
    }
    else { # default to single-ended
        $fastq1 = $fastq;
    }

    my $trimmer;
    switch( lc($trimming_params->{trimmer}) ) {
        case 'cutadapt'    { $trimmer = CoGe::Builder::Trimming::Cutadapt->new($self)     }
        case 'trimgalore'  { $trimmer = CoGe::Builder::Trimming::TrimGalore->new($self)   }
        case 'trimmomatic' { $trimmer = CoGe::Builder::Trimming::Trimmomatic->new($self)  }
        case 'bbduk'       { $trimmer = CoGe::Builder::Trimming::BBDuk->new($self)        }
        default {
            CoGe::Exception::Generic->throw(message => 'Invalid trimmer');
        }
    }

    $trimmer->build(data_files => [$fastq1, $fastq2]);
    $self->add($trimmer);
    push @{$self->fastq}, @{$trimmer->fastq};
}

1;
