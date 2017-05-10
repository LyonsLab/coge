package CoGe::Builder::Load::SRA;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use Hash::Flatten qw(flatten);
use File::Spec::Functions qw(catdir catfile);
use String::ShellQuote;

use CoGe::Accessory::Web qw(url_for);
use CoGe::Accessory::Entrez;
use CoGe::Core::Storage qw(get_sra_cache_path);

use constant MAX_SRA_ACCESSIONS  => 10;
use constant MAX_SRA_EXPERIMENTS => 10;

sub get_name {
    my $self = shift;
    my @accns = map { $_->{path} } @{$self->params->{source_data}};
    return "Import from SRA: " . join(', ', @accns);
}

sub get_site_url {
    my $self = shift;
    return url_for('LoadExperiment.pl', wid => $self->workflow->id);
}

sub build {
    my $self = shift;
    my %opts = @_;
    my $toplevel_accns = $opts{ncbi_accns};
    unless ($toplevel_accns && @$toplevel_accns) {
        CoGe::Exception::Generic->throw(message => 'Missing inputs');
    }

    # Retrieve IDs and metdata from SRA
    my $records = _sra_retrieve($toplevel_accns);
    #warn Dumper $records;

    # Limit the max number of SRA accessions a user can load at one time
    if (@$records >= MAX_SRA_ACCESSIONS && !$self->user->is_poweruser) {
        CoGe::Exception::Generic->throw(message => 'Too many accessions', details => Dumper $toplevel_accns);
    }

    #
    # Build workflow
    #

    foreach my $record (@$records) {
        # Extract metadata from SRA record
        my ($accns, $metadata, $additional_metadata, $read_type) = _extract_metadata($record);
        print STDERR "SRA: Building workflow for '", $metadata->{name}, "'\n";

        # Limit the max number of SRA experiments a user can load at one time
        if (@$accns > MAX_SRA_EXPERIMENTS && !$self->user->is_poweruser) {
            warn "Too many experiments for accessions ", join(', ', @$accns);
            return;
        }

        # Data retrieval is done here instead of Extractor so can be downloaded serially,
        # that is for one experiment at a time (instead of downloading all files up front).
        my @fastq;
        foreach (@$accns) {
            $self->add_to_previous(
                $self->fastq_dump($_, $read_type),
            );
            push @fastq, grep { $_ =~ /\.fastq$/ } @{$self->previous_outputs};
        }

        # Add load experiment pipeline
        $self->params->{metadata} = $metadata;
        $self->params->{additional_metadata} = $additional_metadata;
        $self->params->{read_params}{read_type} = $read_type;
        my $load = CoGe::Builder::Load::Experiment->new($self);
        $load->build(data_files => \@fastq);
        $self->add_to_all($load);
    }
}

sub _sra_retrieve {
    my $accns = shift;

    my @results;
    foreach (@$accns) {
        print STDERR "Querying SRA for $_\n";
        my $id_list = esearch('sra', $_);
        return unless ($id_list && @$id_list);

        print STDERR "   fetch summary for ", join(',', @$id_list), "\n";
        my $records = esummary('sra', $id_list);
        push @results, @$records;
    }

    return \@results;
}

sub _extract_metadata {
    my $record = shift;
    return unless $record;

    # Flatten SRA record metadata
    my $flattened = flatten($record);
    $flattened = {} unless $flattened;
    #warn Dumper $flattened;

    # Extract primary metadata fields
    my $metadata = {
        name       => $flattened->{'Summary.Title'} // '',
        version    => $flattened->{'Experiment.ver'} // 1,
        source_name=> 'NCBI-SRA',
        link       => 'https://www.ncbi.nlm.nih.gov/sra' . ($flattened->{'Experiment.acc'} ? '?term='.$flattened->{'Experiment.acc'} : ''),
        restricted => 0, # data is from a public source
        tags       => [ 'SRA' ]
    };

    # Format other fields as additional metadata
    my @additional_metadata = map +{ type => $_, text => $flattened->{$_} }, keys(%$flattened);
    #warn Dumper \@additional_metadata;

    # Extract accession number(s)
    my @accns;
    if ($record->{Run}->{acc}) { # single run
        push @accns, $record->{Run}->{acc};
    }
    else { # multiple runs
        foreach (keys %{$record->{Run}}) {
            push @accns, $_ if ($_ =~ /^SRR/i);
        }
    }

    # Determine single-ended or paired-end
    my $read_type = 'single';
    if ($record->{'Library_descriptor'}->{'LIBRARY_LAYOUT'}->{'PAIRED'}) { # paired-end
        $read_type = 'paired';
    }

    return (\@accns, $metadata, \@additional_metadata, $read_type);
}

sub fastq_dump {
    my $self = shift;
    my $accn = shift;
    my $read_type = shift // 'single';

    my $dest_path = get_sra_cache_path();

    my $cmd = $self->conf->{FASTQ_DUMP} || 'fastq-dump';
    $cmd .= ' --split-files' if ($read_type eq 'paired');

    my $output_filepath = catfile($dest_path, $accn);

    my (@output_files, @done_files);
    if ($read_type eq 'paired') {
        @output_files = (
            $output_filepath . '_1.fastq',
            $output_filepath . '_2.fastq'
        );
        @done_files = (
            $output_filepath . '_1.fastq.done',
            $output_filepath . '_2.fastq.done'
        );
    }
    else {
        @output_files = ( $output_filepath . '.fastq');
        @done_files   = ( $output_filepath . '.fastq.done' );
    }

    return {
        cmd => "mkdir -p $dest_path && $cmd --outdir $dest_path " . shell_quote($accn) . " && touch " . join(' ', @done_files),
        args => [],
        inputs => [],
        outputs => [
            @output_files,
            @done_files
        ],
        description => "Fetching $accn from NCBI-SRA"
    };
}

__PACKAGE__->meta->make_immutable;

1;
