package CoGe::Builder::Load::Genome;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_upload_path);
use CoGe::Builder::CommonTasks;
use CoGe::Builder::Common::DataRetrieval;

sub get_name {
    my $self = shift;
    my $metadata = $self->params->{metadata};
    my $info;
    $info .= $metadata->{organism} if $metadata->{organism};
    $info .= " (" . $metadata->{name} . ")"  if $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    return "Load Genome \"$info\"";
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $organism_id = $self->params->{organism_id};
    return unless $organism_id;
    my $data = $self->params->{source_data};
    return unless (defined $data && @$data);
    my $metadata = $self->params->{metadata};
    return unless $metadata;
    my $load_id = $self->params->{load_id} || get_unique_id();
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get organism
    my $organism = $self->db->resultset('Organism')->find($organism_id);
    return unless $organism;
    
    #
    # Build workflow
    #
    
    # Create tasks to retrieve files #TODO move to Buildable::pre_build
#    my $upload_dir = get_upload_path($self->user->name, $load_id);
#    my $data_workflow = create_data_retrieval_workflow(upload_dir => $upload_dir, data => $data);
#    $self->add_tasks($data_workflow->{tasks});
#    push @input_files, @{$data_workflow->{outputs}} if ($data_workflow->{outputs});
#    push @ncbi_accns, @{$data_workflow->{ncbi}} if ($data_workflow->{ncbi});
    
    my $request = { #FIXME better way to do this?
        params      => $self->params,
        requester   => $self->requester,
        db          => $self->db,
        user        => $self->user,
        conf        => $self->conf,
        workflow    => $self->workflow,
        staging_dir => $self->staging_dir,
        result_dr   => $self->result_dir,
        outputs     => $self->outputs
    };

    my $dr = CoGe::Builder::Common::DataRetrieval->new($request);
    $dr->build();
    
    # Build steps to add genome
    my @ncbi_accns  = $dr->get_assets('ncbi_accn');
    if (@ncbi_accns) { # NCBI-based load
        $self->add_task(
            $self->load_genome_from_NCBI(
                ncbi_accns => \@ncbi_accns,
                metadata => $metadata,
            )
        );
    }
    else { # File-based load
        my @input_files = $dr->get_assets('data_file');
        
        # Untar/decompress input files
        my (@decompressed, $doJoin);
        foreach my $input_file (@input_files) {
            if ( $input_file =~ /\.tgz|\.tar\.gz$/ ) { # Untar if necessary
                my $task = $self->untar(
                    input_file => $input_file, 
                    output_path => catdir($self->staging_dir, 'untarred')
                );
                $self->add_task($task);
                $input_file = $task->{outputs}->[0] . '/*'; # this is a directory and filespec
                $doJoin = 1;
            }
            elsif ( $input_file =~ /\.gz$/ ) { # Decompress if necessary
                my $task = $self->gunzip( input_file => $input_file );
                $self->add_task($task);
                $input_file = $task->{outputs}->[0];
            }
            push @decompressed, $input_file;
        }
        
        my $concatenated_file;
        if (@decompressed > 1 || $doJoin) {
            # Concatenate all input files into one
            $concatenated_file = catfile($self->staging_dir, 'concatenated_genome.fasta');
            $self->add_task_chain_all(
                $self->join_files(
                    input_files => \@decompressed, # this should work with wildcards from untar task above
                    output_file => $concatenated_file,
#                    done_files => \@done_files
                )
            );
#            push @tasks, $cat_task;
        }
        else {
            $concatenated_file = shift @decompressed; # first and only file
        }
    
        # Sort FASTA by length
        $self->add_task_chain_all(
            $self->sort_fasta( fasta_file => $concatenated_file )
        );
#        push @tasks, $sort_task;
        my $sorted_fasta_file = $self->get_asset('sorted_fasta');
        
        # Validate/process/trim FASTA file
        my $processed_fasta_file = catdir($self->staging_dir, 'genome.faa');
        $self->add_task_chain(
            $self->process_fasta(
                input_file => $sorted_fasta_file,
                output_file => $processed_fasta_file,
            )
        );
#        push @tasks, $process_task;
        
        # Index FASTA file
        $self->add_task_chain(
            create_fasta_index_job(
                fasta => $processed_fasta_file
            )
        );
#        push @tasks, $index_task;
#        my $index_file = $index_task->{outputs}[0];
    
        # Create genome in DB
        $self->add_task_chain(
            $self->load_genome(
                organism_id => $organism->id,
                fasta_file => $processed_fasta_file,
                metadata => $metadata
            )
        );
#        push @tasks, $load_task;
#        push @done_files, $load_task->{outputs}->[1];
    }
    
    return 1;
}

sub sort_fasta {
    my ($self, %params) = @_;
    
    my $fasta_file = $params{fasta_file};
    my $filename = basename($fasta_file);
    my $output_file = "$filename.sorted";
    
    $self->add_asset(sorted_fasta => catfile($self->staging_dir, $output_file));
    
    my $cmd = $self->conf->{SIZESEQ} || 'sizeseq';

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-sequences', $filename, 0],
            ['-descending', 'Y', 0],
            ['-outseq', $output_file, 0]
        ],
        inputs => [
            $fasta_file
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => 'Sorting FASTA file...'
    };
}

sub process_fasta {
    my ($self, %params) = @_;

    my $input_file = $params{input_file};
    my $output_file = $params{output_file};
    
    my $done_file = $input_file . '.processed';
    
    my $cmd = catfile($self->conf->{SCRIPTDIR}, "process_fasta.pl");
    die "ERROR: SCRIPTDIR not specified in config" unless $cmd;

    return {
        cmd => "$cmd -input_fasta_file $input_file -output_fasta_file $output_file && touch $done_file",
        script => undef,
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            $done_file
        ],
        description => "Validating/processing FASTA file...",
    };
}

sub load_genome {
    my ($self, %params) = @_;

    my $metadata    = $params{metadata};
    my $organism_id = $params{organism_id};
    my $fasta_file  = $params{fasta_file};
    
#    my $result_file = get_workflow_results_file($user->name, $wid);

    return {
        cmd => 'perl ' . catfile($self->conf->{SCRIPTDIR}, "load_genome.pl"),
        script => undef,
        args => [
            ['-user_name', $self->user->name, 0],
            ['-wid', $self->workflow->id, 0],
            ['-name', ($metadata->{name} ? shell_quote($metadata->{name}) : '""'), 0],
            ['-desc', ($metadata->{description} ? shell_quote($metadata->{description}) : '""'), 0],
            ['-link', ($metadata->{link} ? shell_quote($metadata->{link}) : '""'), 0],
            ['-version', ($metadata->{version} ? shell_quote($metadata->{version}) :'""'), 0],
            ['-restricted', ($metadata->{restricted} ? 1 : 0), 0],
            ['-source_name', ($metadata->{source_name} ? shell_quote($metadata->{source_name}) : '""'), 0],
            ['-organism_id', $organism_id, 0],
            ['-type_id', ( $metadata->{type} ? shell_quote($metadata->{type}) : 1 ), 0], # default to "unmasked"
            ['-staging_dir', $self->staging_dir, 0],
            ['-fasta_file', shell_quote($fasta_file), 0],
            #['-irods_files', shell_quote($irods_str), 0],
            ['-config', $self->conf->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $fasta_file,
#            @$done_files
        ],
        outputs => [
            catfile($self->staging_dir, "log.done"),
#            $result_file
        ],
        description => "Loading genome ..."
    };
}

sub load_genome_from_NCBI {
    my ($self, %params) = @_;

    my $output_path = catdir($self->staging_dir, "load_genome_from_ncbi");
#    my $result_file = get_workflow_results_file($user->name, $wid);

    my $args = [
        ['-user_name',   $self->user->name,          0],
        ['-wid',         $self->workflow->id,        0],
        ['-staging_dir', "./load_genome_from_ncbi",  0],
        ['-config',      $self->conf->{_CONFIG_PATH}, 0],
        ['-GO',          1,                          0]
    ];
    foreach (@{$params{ncbi_accns}}) {
        push @$args, ['-accn', "'$_'", 0];
    }

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "genbank_genome_loader.pl"),
        script => undef,
        args => $args,
        inputs => [],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
#            $result_file
        ],
        description => "Importing genome from NCBI ..."
    };
}

1;
