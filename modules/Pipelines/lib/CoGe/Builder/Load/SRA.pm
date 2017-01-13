package CoGe::Builder::Load::SRA;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use Hash::Flatten qw(flatten);

use CoGe::Accessory::Web qw(url_for);
use CoGe::Accessory::Entrez;
use CoGe::Core::Storage qw(get_upload_path);
use CoGe::Builder::CommonTasks;

use constant MAX_SRA_ACCESSIONS  => 10;
use constant MAX_SRA_EXPERIMENTS => 10;

sub get_name {
    my $self = shift;
    my @accns = _extract_accns($self->params->{source_data});
    return "Import from SRA: " . join(', ', @accns);
}

sub get_site_url {
    my $self = shift;
    return url_for('LoadExperiment.pl', wid => $self->workflow->id);
}

sub build {
    my $self = shift;
    my %opts = @_;
#    my $input_files = $opts{data_files}; #FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#    # Extract SRA accessions
#    my @accns = _extract_accns($source_data);
#    return unless @accns;
#
#    # Retrieve IDs and metdata from SRA
#    my $records = _sra_retrieve(\@accns);
#
#    # Limit the max number of SRA accessions a user can load at one time
#    if (@$records >= MAX_SRA_ACCESSIONS && !$self->user->is_poweruser) {
#        CoGe::Exception::Generic->throw(message => 'Too many accessions', details => Dumper \@accns);
#    }
#
#    #
#    # Build workflow
#    #
#
#    my $wait_file;
#    foreach my $record (@$records) {
#        # Format metadata
#        my ($data, $metadata, $additional_metadata, $read_type) = _convert_metadata($record);
#        print STDERR 'SRA: Building workflow for "', $metadata->{name}, '"\n';
#
#        # Limit the max number of SRA experiments a user can load at one time
#        if (@$data > MAX_SRA_EXPERIMENTS && !$self->user->is_poweruser) {
#            warn "Too many experiments for accessions ", join(', ', @accns);
#            return;
#        }
#
#        # Add load experiment pipeline
#        $self->params->{metadata} = $metadata;
#        $self->params->{additional_metadata} = $additional_metadata;
#        $self->params->{read_params}{read_type} = $read_type;
#        $self->params->{source_data} = $data;
#
#        my $exp = CoGe::Builder::Load::Experiment->new();
#        $exp->build();
#
#        # Add wait task to serialize experiment loads -- #FIXME kludge until this can be generalized and incorporated into Buildable
#        $exp->add_task_chain_all(
#            $self->create_wait()
#        );
##        my @all_outputs = $self->workflow->get_outputs();
##        my $wait_task = $self->create_wait();
##        push @{$wait_task->{inputs}}, @all_outputs;
##        $self->workflow->add_job($wait_task);
##        $wait_file = $wait_task->{outputs}->[0];
#    }
}

sub _extract_accns {
    my $source_data = shift;

    # Extract SRA accessions and verify only SRA data types present
    my @accns;
    foreach (@$source_data) {
        if (!$_->{type} || lc($_->{type}) ne 'sra') {
            warn "Data type is undefined or not 'sra'";
            return;
        }

        if (!$_->{path}) {
            warn "Data path is undefined";
            return;
        }

        push @accns, $_->{path};
    }

    return wantarray ? @accns : \@accns;
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

sub _convert_metadata { #TODO rename
    my $record = shift;
    return unless $record;

    # Flatten SRA record metadata
    my $flattened = flatten($record);
    $flattened = {} unless $flattened;
    #warn Dumper $record;
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
    my $read_type;
    if ($record->{'Library_descriptor'}->{'LIBRARY_LAYOUT'}->{'PAIRED'}) { # paired-end
        $read_type = 'paired';
    }
    # Else default to single-ended

    # Create data element
    my @data;
    foreach (@accns) {
        push @data, { file_type => 'sra', type => 'sra', path => $_ };
    }

    return (\@data, $metadata, \@additional_metadata, $read_type);
}

__PACKAGE__->meta->make_immutable;

1;
