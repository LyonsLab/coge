package CoGe::Builder::Tools::Pseudoassembly;

use Moose;
extends 'CoGe::Builder::Buildable';

use File::Spec::Functions;
use Digest::MD5 qw(md5_hex);
use CoGe::Accessory::Web qw(download_url_for);

sub get_name {
	my $self = shift;
    return 'Generate pseudo assembly for ' . $self->request->genome1->info . ' v. ' . $self->request->genome2->info;
}

sub build {
	my $self = shift;

	my $input = $self->params->{input};
	my $flip  = $self->params->{flip} // 0;

	my $gid1  = $self->request->genome1->id;
	my $gid2  = $self->request->genome2->id;

	# Swap order if gid2 < gid1
	( $gid1, $gid2 ) = ( $gid2, $gid1 ) if $gid2 lt $gid1;

	my $filename = qq($gid1-$gid2-) . md5_hex($input) . '.' . $flip . '.tar.gz';
	my $output = catfile($self->conf->{DIAGSDIR}, $gid1, $gid2, 'assembly', $filename);

	$self->add({
		cmd  => catfile( $self->conf->{SCRIPTDIR}, 'synmap/order_contigs_to_chromosome.pl' ),
		args => [
			[ "-cfg",    $self->conf->{_CONFIG_PATH}, 1 ],
			[ "-input",  $input,                  1 ],
			[ "-output", $output,                 1 ],
			[ "-flip",   $flip,                   1 ],
		],
		inputs  => [ $input, $self->conf->{_CONFIG_PATH} ],
		outputs => [$output],
		description => "Generating pseudo assembly"
	});

	my $dir = $self->conf->{COGEDIR};
	my $url = $self->conf->{URL};
	$output =~ s/$dir/$url/; #TODO move this to API download endpoint

	$self->add_to_previous(
		$self->add_result(
			result   => {
				type => 'url',
				name => 'pseudoassembly',
				path => $output
			}
		)
	);
}

__PACKAGE__->meta->make_immutable;

1;