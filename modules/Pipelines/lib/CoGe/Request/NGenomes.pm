package CoGe::Request::NGenomes;

use Moose;
extends 'CoGe::Request::Request';

has genomes => (is => 'rw', isa => 'ArrayRef', default => sub{ [] });

sub is_valid {
    my $self = shift;

	for (my $i = 1; my $id = $self->parameters->{'genome_id' . $i}; $i++) {
        my $genome = $self->db->resultset("Genome")->find($id);
		return 0 unless $genome;

        push @{$self->genomes}, $genome;
	}

    return scalar(@{$self->genomes}) > 2 ? 1 : 0;
}

sub has_access {
    my $self = shift;

    foreach (@{$self->genomes}) {
        return 0 if ($_->restricted && (!$self->user || !$self->user->has_access_to_genome($_)));
    }

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
