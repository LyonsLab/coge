package CoGe::Builder::Tools::NucCounter;

use Moose;
with qw(CoGe::Builder::Buildable);

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use File::Spec::Functions;

sub pre_build { # override superclass method for reusable workflow ID, custom site_url, and custom workflow paths
	my ($self, %params) = @_;

	# Initialize workflow -- NOTE: init => 0 means that a previous identical workflow will be reused when submitted
    $self->workflow( $params{jex}->create_workflow(name => $self->get_name, init => 0 ) );
    return unless $self->workflow;
}

sub build {
	my $self = shift;

    my $gid = $self->params->{'gid'};
    my $chr = $self->params->{'chr'};
    my $dir = catfile($self->conf->{SECTEMPDIR}, "downloads/genome", $gid);

    my $fasta = catfile($dir, $gid . '_' . $chr . '.faa');
    $self->add_task({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'generate_chr_fasta.pl'),
        args        => [[ 'gid', $gid, 0 ], [ 'chr', $chr, 0 ]],
        outputs     => [$fasta],
        description => "Generating chromosome sequence",        
    });

    my $filename = $gid . '_' . $chr . '_out.txt';
    my $output = catfile($dir, $filename);
    $self->add_task({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'nuccounter.py') . ' ' . $fasta,
        inputs      => [$fasta],
        outputs     => [$output],
        description => "Generating nucleotide sliding window percentages",        
    });

    if ($self->params->{'irods'}) {
        my $irods_base = irods_get_base_path($self->user->name);
        $self->add_task(
            $self->export_to_irods(
                src_file  => $output,
                dest_file => catfile($irods_base, $filename)
            )
        );
    }

	return 1;
}

sub get_name {
	my $self = shift;
    return 'NucCounter | ' . $self->params->{'gid'} . ' | ' . $self->params->{'chr'};
}