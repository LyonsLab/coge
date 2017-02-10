package CoGe::Builder::Trimming::BBDuk;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my %opts = @_;
    my ($fastq1, $fastq2) =  @{$opts{data_files}}; # fastq2 is undef for single-ended
    unless ($fastq1 && @$fastq1) {
        CoGe::Exception::MissingField->throw(message => 'Missing fastq');
    }

    my $genome = $self->request->genome;

    # Reheader the fasta file
    $self->add(
        $self->reheader_fasta($genome->id)
    );
    my $reheader_fasta = $self->previous_output;

    # Index the fasta file
    $self->add_to_previous(
        $self->index_fasta($self->previous_output)
    );

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    if ($read_type eq 'single') { # single-ended
        # Create BBDuk task for each file
        foreach (@$fastq1) {
            $self->add(
                $self->bbduk($_, $reheader_fasta)
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create BBDuk task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            $self->add(
                $self->bbduk([ $fastq1->[$i], $fastq2->[$i] ], $reheader_fasta)
            );
            push @{$self->fastq}, grep { /\.fastq$/ } @{$self->previous_outputs};
        }
    }
}

sub bbduk {
    my $self = shift;
    my $fastq = shift; # single fastq file or two paired-end fastq files
    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');
    my $fasta = shift;

    my $trimming_params = $self->params->{trimming_params} // {};
    my $read_params = $self->params->{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    my @outputs = map { catfile($self->staging_dir, to_filename($_) . '.trimmed.fastq') } @$fastq;

    unless ($self->conf->{BBMAP}) {
        CoGe::Exception::Generic->throw(message => 'Missing BBMAP in configuration file');
    }

    my $cmd = join(' ',
        'nice',
        catfile($self->conf->{BBMAP}, 'bbduk.sh'),
        qq[in=$fastq->[0]],
        qq[out=$outputs[0]],
        qq[ref=$fasta]
    );

    if ($read_type eq 'paired') { # paired-end
        $cmd .= join(' ',
            qq[in2=$fastq->[1]],
            qq[out2=$outputs[1]]
        );
    }

    return {
        cmd         => $cmd,
        args        => [],
        inputs      => [ @$fastq ],
        outputs     => \@outputs,
        description => 'Trimming (BBDuk) '.join(', ', map { basename($_) } @$fastq)
    };
}

1;