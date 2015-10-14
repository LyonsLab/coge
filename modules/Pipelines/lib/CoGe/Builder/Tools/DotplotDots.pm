package CoGe::Builder::Tools::DotplotDots;

use Moose;

use CoGe::Core::Storage qw( get_workflow_paths );

sub build {
	my $self = shift;

	# Validate inputs
	my $genome_id1 = $self->params->{genome_id1};
	return unless $genome_id1;
	my $genome_id2 = $self->params->{genome_id2};
	return unless $genome_id2;
	my $ksfile = $self->params->{ksfile};
	return unless $ksfile;

	$self->workflow( $self->jex->create_workflow( name => "DotplotDots", init => 1 ) );
	return unless ( $self->workflow && $self->workflow->id );
	my ( $staging_dir, $result_dir ) = get_workflow_paths( $self->user->name, $self->workflow->id );
	$self->workflow->logfile( catfile( $result_dir, "debug.log" ) );

	my $DIAGSDIR      = $self->conf->{DIAGSDIR};
	my $SCRIPTDIR     = $self->conf->{SCRIPTDIR};
	my ( $dir1, $dir2 ) = sort ( $genome_id1, $genome_id2 );

	$self->workflow->add_job(
		{
			cmd         => catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile,
			outputs     => [catfile($DIAGSDIR, $dir1, $dir2, 'log.json')],
			description => "running dotplot_dots...",
		}
	);

	return 1;
}

with qw(CoGe::Builder::Buildable);

1;
