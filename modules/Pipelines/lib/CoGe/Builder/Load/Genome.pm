package CoGe::Builder::Load::Genome;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path url_for);
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;
use CoGe::Exception::ItemNotFound;

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
    my %opts = @_;
    my $input_files = $opts{data_files};
    my $input_dir   = $opts{data_dir};
    my $ncbi_accns  = $opts{ncbi_accns};
    unless (@$input_files || @$input_dir || @$ncbi_accns) {
        CoGe::Exception::MissingField->throw(message => "Missing inputs");
    }
    
    # Validate inputs
    my $organism_id = $self->params->{organism_id};
    unless ($organism_id) {
        CoGe::Exception::MissingField->throw(message => "Missing organism_id");
    }
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }

    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get organism
    my $organism = $self->db->resultset('Organism')->find($organism_id);
    unless ($organism) {
        CoGe::Exception::ItemNotFound->throw(type => 'organism', id => $organism_id);
    }
    
    #
    # Build workflow
    #
    
    # Build steps to add genome
    if (@$ncbi_accns) { # NCBI-based load
        $self->add(
            $self->load_genome_from_NCBI(
                ncbi_accns => $ncbi_accns,
                metadata => $metadata,
            )
        );
    }
    else { # File-based load
        my ($fasta_file) = @$input_files; # first file (in case just one);
        if (@$input_files > 1 || $input_dir) { # multiple FASTA files
            # Concatenate all input files into one
            $self->add_to_all(
                $self->join_files(
                    input_files => $input_files,
                    input_dir => $input_dir,
                    output_file => catfile($self->staging_dir, 'concatenated_genome.fasta'),
                )
            );
            $fasta_file = $self->previous_output;
        }
    
        # Sort FASTA by length (for trimming in next step)
        $self->add_to_all(
            $self->sort_fasta( fasta_file => $fasta_file )
        );
        
        # Validate/process/trim FASTA file
        my $processed_fasta_file = catdir($self->staging_dir, 'genome.faa');
        $self->add_to_previous(
            $self->process_fasta(
                input_file => $self->previous_output,
                output_file => $processed_fasta_file,
            )
        );
        
        # Index processed FASTA file
        $self->add_to_previous(
            $self->index_fasta($self->previous_output)
        );
    
        # Create genome in DB
        $self->add_to_previous(
            $self->load_genome(
                organism_id => $organism->id,
                fasta_file => $processed_fasta_file,
            )
        );
    }
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
    my $organism_id = $params{organism_id};
    my $fasta_file  = $params{fasta_file};

    my $metadata    = $self->params->{metadata};
    
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

__PACKAGE__->meta->make_immutable;

1;
