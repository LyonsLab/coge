package CoGe::Builder::Load::Genome;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Accessory::Web qw(url_for);
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

sub get_site_url {
    my $self = shift;
    return url_for('LoadGenome.pl', wid => $self->workflow->id);
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

    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get organism
    my $organism = $self->db->resultset('Organism')->find($organism_id);
    return unless $organism;
    
    #
    # Build workflow
    #
    
    # Create tasks to retrieve files #TODO move to pre_build()
    my $dr = CoGe::Builder::Common::DataRetrieval->new({ #FIXME better way to pass these args? See Moose constructors
#        params      => $self->params,
#        db          => $self->db,
#        user        => $self->user,
#        conf        => $self->conf,
        request     => $self->request,
        workflow    => $self->workflow,
        staging_dir => $self->staging_dir,
        result_dr   => $self->result_dir,
        outputs     => $self->outputs
    });
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
        my $input_dir   = $dr->get_asset('data_dir');
        
        my ($fasta_file) = @input_files; # first file (in case just one);
        if (@input_files > 1 || $input_dir) {
            # Concatenate all input files into one
            $self->add_task_chain_all(
                $self->join_files(
                    input_files => \@input_files, 
                    input_dir => $input_dir,
                    output_file => catfile($self->staging_dir, 'concatenated_genome.fasta'),
                )
            );
            $fasta_file = $self->previous_output();
        }
    
        # Sort FASTA by length (for trimming in next step)
        $self->add_task_chain_all(
            $self->sort_fasta( fasta_file => $fasta_file )
        );
        
        # Validate/process/trim FASTA file
        my $processed_fasta_file = catdir($self->staging_dir, 'genome.faa');
        $self->add_task_chain(
            $self->process_fasta(
                input_file => $self->previous_output(),
                output_file => $processed_fasta_file,
            )
        );
        
        # Index processed FASTA file
        $self->add_task_chain(
            create_fasta_index_job(
                fasta => $self->previous_output()
            )
        );
    
        # Create genome in DB
        $self->add_task_chain(
            $self->load_genome(
                organism_id => $organism->id,
                fasta_file => $processed_fasta_file,
                metadata => $metadata
            )
        );
    }
    
    return 1;
}

# mdb removed 10/12/16 -- replaced below for COGE-721
#sub sort_fasta {
#    my ($self, %params) = @_;
#
#    my $fasta_file = $params{fasta_file};
#    my $filename = basename($fasta_file);
#    my $output_file = "$filename.sorted";
#
#    my $cmd = $self->conf->{SIZESEQ} || 'sizeseq';
#
#    return {
#        cmd => $cmd,
#        script => undef,
#        args => [
#            ['-sequences', $filename, 0],
#            ['-descending', 'Y', 0],
#            ['-outseq', $output_file, 0]
#        ],
#        inputs => [
#            $fasta_file
#        ],
#        outputs => [
#            catfile($self->staging_dir, $output_file)
#        ],
#        description => 'Sorting FASTA file'
#    };
#}
sub sort_fasta {
    my ($self, %params) = @_;

    my $fasta_file = $params{fasta_file};
    my $filename = basename($fasta_file);
    my $output_file = "$filename.sorted";

    my $cmd = catfile($self->conf->{SCRIPTDIR}, 'sort_fasta.pl');

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['',  $fasta_file,  1],
            ['>', $output_file, 0]
        ],
        inputs => [
            $fasta_file
        ],
        outputs => [
            catfile($self->staging_dir, $output_file)
        ],
        description => 'Sorting FASTA file'
    };
}

sub process_fasta {
    my ($self, %params) = @_;

    my $input_file = $params{input_file};
    my $output_file = $params{output_file};
    
    my $done_file = $input_file . '.processed';
    
    my $cmd = catfile($self->conf->{SCRIPTDIR}, "process_fasta.pl");

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
        description => "Validating/processing FASTA file",
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
        description => "Loading genome"
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
        description => "Importing genome from NCBI"
    };
}

1;
