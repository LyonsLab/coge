package CoGe::Builder::Load::BatchExperiment;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub get_name {
    my $self = shift;
    my $metadata = $self->params->{metadata};
    my $info = '"' . $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    return "Load batch " . $info;
}

sub build {
    my $self = shift;
    my %opts = @_;
    my $input_files = $opts{data_files};

    my $genome = $self->request->genome;
    
    # Validate inputs
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }
    my $load_id = $self->params->{load_id} || get_unique_id();
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Determine file type if not set
    my $data = $self->params->{source_data};
    my $file_type = $data->[0]->{file_type}; # type of first data file
    ($file_type) = detect_data_type($file_type, $data->[0]->{path}) unless $file_type;
    
    #
    # Build workflow
    #

    # Add load batch task
    $self->add(
        $self->load_batch(
            gid => $genome->id,
            nid => $self->params->{notebook_id},
            input_files => $input_files
        )
    );
}

sub load_batch {
    my $self = shift;
    my %opts = @_;
    my $nid = $opts{nid};
    my $gid = $opts{gid};
    my $files = $opts{input_files};

    my $file_str = join(',', @$files);

    my $metadata = $self->params->{metadata};

    my $args = [
        ['-user_name',   $self->user->name, 0],
        ['-name',        shell_quote($metadata->{name}), 0],
        ['-desc',        shell_quote($metadata->{description}), 0],
        ['-gid',         $gid, 0],
        ['-wid',         $self->workflow->id, 0],
        ['-staging_dir', $self->staging_dir, 0],
        ['-files',       "'".$file_str."'", 0],
        ['-config',      $self->conf->{_CONFIG_PATH}, 0]
    ];
    push $args, ['-nid', $nid, 0] if ($nid);

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "load_batch.pl"),
        args => $args,
        inputs => [
            @$files
        ],
        outputs => [
            [$self->staging_dir, 1],
            catdir($self->staging_dir, 'log.done')
        ],
        description => "Loading batch experiments"
    };
}

__PACKAGE__->meta->make_immutable;

1;
