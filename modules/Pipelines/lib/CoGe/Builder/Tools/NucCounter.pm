package CoGe::Builder::Tools::NucCounter;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use File::Spec::Functions;

sub build {
	my $self = shift;

    my $gid = $self->params->{'gid'};
    my $chr = $self->params->{'chr'};
    my $dir = catfile($self->conf->{SECTEMPDIR}, "downloads/genome", $gid);

    my $fasta = catfile($dir, $gid . '_' . $chr . '.faa');
    $self->add({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'generate_chr_fasta.pl'),
        args        => [[ 'gid', $gid, 0 ], [ 'chr', $chr, 0 ]],
        outputs     => [$fasta],
        description => "Generating chromosome sequence",        
    });

    my $filename = $gid . '_' . $chr . '_out.txt';
    my $output = catfile($dir, $filename);
    $self->add({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'nuccounter.py') . ' ' . $fasta,
        inputs      => [$fasta],
        outputs     => [$output],
        description => "Generating nucleotide sliding window percentages",        
    });

    if ($self->params->{'irods'}) {
        my $irods_base = irods_get_base_path($self->user->name);
        $self->add(
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

1;