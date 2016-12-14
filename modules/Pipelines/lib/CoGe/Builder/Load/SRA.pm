package CoGe::Builder::Load::SRA;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename);
use Hash::Flatten qw(flatten);

use CoGe::Accessory::Web qw(url_for);
use CoGe::Accessory::Entrez;
use CoGe::Core::Storage qw(get_upload_path);
use CoGe::Builder::CommonTasks;

use constant MAX_SRA_ACCESSIONS  => 10;
use constant MAX_SRA_EXPERIMENTS => 10;

sub get_name {
    my $self = shift;
    return "Import from SRA"; #TODO add accn(s)
}

sub get_site_url {
    my $self = shift;
    return url_for('LoadExperiment.pl', wid => $self->workflow->id);
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $gid = $self->params->{genome_id} || $self->params->{gid};
    return unless $gid;
    my $source_data = $self->params->{source_data};
    return unless (defined $source_data && @$source_data);
#    my $metadata = $self->params->{metadata}; # primary metadata is exctracted from SRA record
#    return unless $metadata;
    my $additional_metadata = $self->params->{additional_metadata}; # optional

    # Convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
#    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get genome
    my $genome = $self->db->resultset('Genome')->find($gid);
    return unless $genome;

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

    # Retrieve IDs and metdata from SRA
    my $records = sra_retrieve(\@accns);

    # Limit the max number of SRA accessions a user can load at one time
    if (@$records >= MAX_SRA_ACCESSIONS && !$self->user->is_poweruser) {
        warn "Too many accessions ", join(', ', @accns);
        return;
    }

    #
    # Build workflow
    #

    my $wait_file;
    foreach my $record (@$records) {
        # Format metadata
        my ($data, $metadata, $additional_metadata, $read_params) = convert_metadata($record);
        print STDERR 'SRA: Building workflow for "', $metadata->{name}, '"\n';

        # Limit the max number of SRA experiments a user can load at one time
        if (@$data > MAX_SRA_EXPERIMENTS && !$self->user->is_poweruser) {
            warn "Too many experiments for accessions ", join(', ', @accns);
            return;
        }

        my $expBuilder = CoGe::Builder::Load::Experiment->new({
            params      => {
                genome_id           => $gid,
                metadata            => $metadata,
                additional_metadata => $additional_metadata,
                source_data         => $data,
                read_params         => $read_params
            },
            db          => $self->db,
            user        => $self->user,
            conf        => $self->conf,
            workflow    => $self->workflow,
            staging_dir => $self->staging_dir,
            result_dr   => $self->result_dir,
            inputs      => ($wait_file ? [ $wait_file ] : undef),
            #outputs     => $self->outputs
        });
        $expBuilder->build();

        # Add wait task to serialize experiment loads -- #FIXME kludge until this can be generalized and incorporated into Buildable
#        $expBuilder->add_task_chain_all( # won't work because Experiment.pm doesn't use Buildable::add_job
#            $self->create_wait()
#        );
#        $wait_file = $self->previous_output;
        my @all_outputs = $self->workflow->get_outputs();
        my $wait_task = $self->create_wait();
        push @{$wait_task->{inputs}}, @all_outputs;
        $self->workflow->add_job($wait_task);
        $wait_file = $wait_task->{outputs}->[0];

        last; #debug
    }
    
    return 1;
}

sub sra_retrieve {
    my $accns = shift;

    my @results;
    foreach (@$accns) {
        print STDERR "Querying SRA for $_\n";
        my $id_list = esearch('sra', $_);
        return unless ($id_list && @$id_list);

        foreach (@$id_list) { #FIXME speed this up by fetching all in a single query using command-separated list of IDs
            print STDERR "   fetch summary for $_\n";
            my $record = esummary('sra', $_);
            push @results, $record;
        }
    }

    return \@results;
}

sub convert_metadata { #TODO rename
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
    warn Dumper \@additional_metadata;

    # Extract accession number(s) and determine single-ended or paired-end
    my (@accns, $read_params);
    if ($record->{Run}->{acc}) { # single run
        push @accns, $record->{Run}->{acc};
    }
    else { # multiple runs
        foreach (keys %{$record->{Run}}) {
            push @accns, $_ if ($_ =~ /^SRR/i);
        }

        if ($record->{'Library_descriptor'}->{'LIBRARY_LAYOUT'}->{'PAIRED'}) { # paired-end
            $read_params->{read_type} = 'paired';
        }
        # Else default to single-ended
    }

    # Create data element
    my @data;
    foreach (@accns) {
        push @data, { file_type => 'sra', type => 'sra', path => $_ };
    }

    return (\@data, $metadata, \@additional_metadata, $read_params);
}

1;
