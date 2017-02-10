package CoGe::Request::NGenomes;

use Moose;
extends 'CoGe::Request::Request';

use CoGe::Exception::Generic;
use CoGe::Exception::AccessDenied;
use CoGe::Exception::ItemNotFound;

has genomes => (is => 'rw', isa => 'ArrayRef', default => sub{ [] });

sub is_valid {
    my $self = shift;

	for (my $i = 1; my $id = $self->parameters->{'genome_id' . $i}; $i++) {
        my $genome = $self->db->resultset("Genome")->find($id);
        unless ($genome) {
            CoGe::Exception::ItemNotFound->throw(type => 'genome', id => $id);
        }

        push @{$self->genomes}, $genome;
	}

    if (scalar(@{$self->genomes}) < 2) {
        CoGe::Exception::Generic->throw(message => 'At least 2 genomes must be specified');
    }

    return 1;
}

sub has_access {
    my $self = shift;

    foreach my $genome (@{$self->genomes}) {
        if ($genome->restricted && (!$self->user || !$self->user->has_access_to_genome($genome))) {
            CoGe::Exception::AccessDenied->throw(type => 'genome', id => $genome->id);
        }
    }

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
