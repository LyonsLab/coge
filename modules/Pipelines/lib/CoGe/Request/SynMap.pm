package CoGe::Request::SynMap;

use Moose;
extends 'CoGe::Request::Request';

sub is_valid {
    my $self = shift;
    my $id;
    my $i;
	for ($i=1; $id=$self->parameters->{'genome_id' . $i}; $i++) {
		return 0 unless $self->db->resultset("Genome")->find($id);
	}
    return $i > 2 ? 1 : 0;
}

sub has_access {
    my $self = shift;
    my $id;
	for (my $i=1; $id=$self->parameters->{'genome_id' . $i}; $i++) {
	    my $genome = $self->db->resultset("Genome")->find($id);
	    return 0 if ($genome->restricted && (!$self->user || !$self->user->has_access_to_genome($genome)));
	}
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
