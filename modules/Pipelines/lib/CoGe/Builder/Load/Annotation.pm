package CoGe::Builder::Load::Annotation;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub get_name {
    my $self = shift;
    my $genome = $self->request->genome;
    my $metadata = $self->params->{metadata};

    my $info = '"';
    $info .= $genome->organism->name;
    $info .= " (" . $metadata->{name} . ")"  if $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    return "Load Annotation " . $info;
}

sub build {
    my $self = shift;
    my %opts = @_;
    my ($gff_file) = @{$opts{data_files}};
    unless ($gff_file) {
        CoGe::Exception::Generic->throw(message => 'Missing gff');
    }

    # Validate inputs
    my $genome = $self->request->genome;
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }

    #
    # Build workflow
    #

    $self->add(
        $self->load_annotation(
            gid => $genome->id,
            input_file => $gff_file
        )
    );
}

sub load_annotation {
    my $self = shift;
    my %opts = @_;
    my $gid = $opts{gid};
    my $input_file = $opts{input_file};

    my $output_path = catdir($self->staging_dir, "load_annotation");
#    my $result_file = get_workflow_results_file($user->name, $wid);

    my $metadata = $self->params->{metadata};

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "load_annotation.pl"),
        args => [
            ['-user_name',   $self->user->name,   0],
            ['-wid',         $self->workflow->id, 0],
            ['-name',        ($metadata->{name} ?        shell_quote($metadata->{name})    : qq['']),     0],
            ['-desc',        ($metadata->{description} ? shell_quote($metadata->{description}) : qq['']), 0],
            ['-link',        ($metadata->{link} ?        shell_quote($metadata->{link})    : qq['']),     0],
            ['-version',     ($metadata->{version} ?     shell_quote($metadata->{version}) : qq['']),     0],
            ['-source_name', ($metadata->{source} ?      shell_quote($metadata->{source})  : qq['']),     0],
            ['-gid',         $gid, 0],
            ['-staging_dir', "./load_annotation",         0],
            ['-data_file',   shell_quote($input_file),    0],
            ['-config',      $self->conf->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $input_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
#            $result_file
        ],
        description => "Loading annotation"
    };
}


__PACKAGE__->meta->make_immutable;

1;
