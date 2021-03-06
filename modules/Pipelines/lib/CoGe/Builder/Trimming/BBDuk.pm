package CoGe::Builder::Trimming::BBDuk;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

# Settings
my $JAVA_MAX_MEM = '20g'; # Java max memory size in GB -- value set for TERRA data

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
                $self->bbduk($_, $reheader_fasta),
                qq[$_.done] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create BBDuk task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my ($f1, $f2) = ($fastq1->[$i], $fastq2->[$i]);
            $self->add(
                $self->bbduk([ $fastq1->[$i], $fastq2->[$i] ], $reheader_fasta),
                [ qq{$f1.done}, qq{$f2.done} ] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, @{$self->previous_outputs};
        }
    }
}

sub bbduk {
    my $self = shift;
    my $fastq = shift; # single fastq file or two paired-end fastq files
    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');
    my $fasta = shift;

    my $trimming_params = $self->params->{trimming_params} // {};
#    my $k         = $trimming_params->{k} // 0;
#    my $mink      = $trimming_params->{mink} // -1;
#    my $hdist     = $trimming_params->{hdist} // 0;
#    my $tpe       = $trimming_params->{tpe} // 'f';
#    my $tbo       = $trimming_params->{tbo} // 'f';
#    my $qtrim     = $trimming_params->{qtrim} // 'f';
#    my $trimq     = $trimming_params->{trimq} // 6;
#    my $minlength = $trimming_params->{minlength} // 10;
    my $adapters  = $trimming_params->{adapters} // 'none';

    my $read_params = $self->params->{read_params} // {};
#    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    my @outputs = map { catfile($self->staging_dir, basename(remove_fastq_ext($_) . '.trimmed.fastq' . to_compressed_ext($_) . '.done')) } @$fastq;

    unless ($self->conf->{BBMAP}) {
        CoGe::Exception::Generic->throw(message => 'Missing BBMAP in configuration file');
    }

    # Set inputs/outputs
    my $cmd = join(' ',
        'nice',
        catfile($self->conf->{BBMAP}, 'bbduk2.sh'),
        qq[-Xmx$JAVA_MAX_MEM],
        qq[in=$fastq->[0]],
        qq[out=$outputs[0]],
        qq[ref=$fasta]
    );

    if ($read_type eq 'paired') { # paired-end
        $cmd .= ' ' . join(' ',
            qq[in2=$fastq->[1]],
            qq[out2=$outputs[1]]
        );
    }

    # Set threads
    $cmd .= ' threads=' . $self->NUM_CPUS;

    # Set user parameters
    $cmd .= ' ' . join(' ',
        map { qq[$_=$trimming_params->{$_}] } ('k', 'mink', 'hdist', 'tpe', 'tbo', 'qtrim', 'trimq', 'minlength')
    );

    if ($adapters && $adapters ne 'none') {
        my $adapter_file_path = catfile($self->conf->{RESOURCEDIR}, 'adapters.fa');
        if ($adapters eq 'r' || $adapters eq 'both') {
            $cmd .= ' rref=' . $adapter_file_path;
        }
        if ($adapters eq 'l' || $adapters eq 'both') {
            $cmd .= ' lref=' . $adapter_file_path;
        }
    }

    $cmd .= ';touch ' . qq[$outputs[0]];
    $cmd .= ';touch ' . qq[$outputs[1]] if $read_type eq 'paired';

    return {
        cmd         => $cmd,
        args        => [],
        inputs      => [
            @$fastq,
            $fasta
        ],
        outputs     => \@outputs,
        description => 'Trimming (BBDuk) '. fastq_description($fastq, $read_type)
    };
}

1;