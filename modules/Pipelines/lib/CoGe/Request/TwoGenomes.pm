package CoGe::Request::TwoGenomes;

use Moose;
extends 'CoGe::Request::Request';

has genome1 => (is => 'rw', isa => 'CoGeX::Result::Genome');
has genome2 => (is => 'rw', isa => 'CoGeX::Result::Genome');

sub is_valid {
    my $self = shift;

    my $genome_id1 = $self->parameters->{genome_id1};
    return unless $genome_id1;
    my $genome_id2 = $self->parameters->{genome_id2};
    return unless $genome_id2;

    my $genome1 = $self->db->resultset("Genome")->find($genome_id1);
    return unless defined $genome1;
    my $genome2 = $self->db->resultset("Genome")->find($genome_id2);
    return unless defined $genome2;

    $self->genome1($genome1);
    $self->genome2($genome2);

    return 1;
}

sub has_access {
    my $self = shift;
    return unless defined $self->{user};

    unless ($self->user->has_access_to_genome($self->genome1)) {
        CoGe::Exception::AccessDenied->throw(type => 'genome', id => $self->genome1->id);
    }
    unless ($self->user->has_access_to_genome($self->genome2)) {
        CoGe::Exception::AccessDenied->throw(type => 'genome', id => $self->genome2->id);
    }

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
