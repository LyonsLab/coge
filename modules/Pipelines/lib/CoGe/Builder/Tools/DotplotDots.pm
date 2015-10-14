package CoGe::Builder::Tools::DotplotDots;

use Moose;

use File::Spec::Functions;

sub build {
	my $self = shift;

	# Validate inputs
	my $genome_id1 = $self->params->{genome_id1};
	return unless $genome_id1;
	my $genome_id2 = $self->params->{genome_id2};
	return unless $genome_id2;
	my $ksfile = $self->params->{ksfile};
	return unless $ksfile;

	my $DIAGSDIR      = $self->conf->{DIAGSDIR};
	my $SCRIPTDIR     = $self->conf->{SCRIPTDIR};
	my ( $dir1, $dir2 ) = sort ( $genome_id1, $genome_id2 );

	$self->workflow->add_job({
		cmd         => catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile,
		outputs     => [catfile($DIAGSDIR, $dir1, $dir2, 'log.json')],
		description => "running dotplot_dots...",
	});

	return 1;
}

sub get_name {
	return 'DotplotDots';
}

with qw(CoGe::Builder::Buildable);

1;
