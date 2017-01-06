package CoGe::Builder::Tools::NucCounter;

use Moose;
with qw(CoGe::Builder::Buildable);

sub pre_build { # override superclass method for reusable workflow ID, custom site_url, and custom workflow paths
	my ($self, %params) = @_;

	# Initialize workflow -- NOTE: init => 0 means that a previous identical workflow will be reused when submitted
    $self->workflow( $params{jex}->create_workflow(name => $self->get_name, init => 0 ) );
    return unless $self->workflow;
}

sub build {
	my $self = shift;
    my $workflow = $self->workflow;

    my $dir = catfile($config->{SECTEMPDIR}, "downloads/genome", $gid);

    my $fasta = catfile($dir, $gid . '_' . $chr . '.faa');
    $workflow->add_task({
        script      => catfile($self->conf->{SCRIPTDIR}, 'generate_chr_fasta.pl'),
        args        => [[ 'gid', $gid, 0 ], [ 'chr', $chr, 0 ]],
        outputs     => [$fasta],
        description => "Fetching chromosome sequence",        
    });

     my $output = catfile($dir, $gid . '_' . $chr . '_out.txt');
    $workflow->add_task({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'nuccounter.py') . ' ' . $fasta,
        inputs      => [$fasta],
        outputs     => [$fasta1],
        description => "Generating nucleotide sliding window percentages",        
    });

	return 1;
}