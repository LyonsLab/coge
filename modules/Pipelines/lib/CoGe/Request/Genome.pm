package CoGe::Request::Genome;

use Moose;
extends 'CoGe::Request::Request';

use CoGe::Exception::ItemNotFound

has genome => (is => 'rw', isa => 'CoGeX::Result::Genome');

sub is_valid { # called first
    my $self = shift;

    my $gid = $self->parameters->{gid} || $self->parameters->{genome_id};
    unless ($gid) {
        CoGe::Exception::MissingField->throw(message => "Missing genome_id");
    }

    my $genome = $self->db->resultset("Genome")->find($gid);
    unless ($genome) {
        CoGe::Exception::ItemNotFound->throw(type => 'genome', id => $gid);
    }
    $self->genome($genome);

    return 1;
}

sub has_access { # called second
    my $self = shift;

    if ($self->genome->restricted && (!$self->user || !$self->user->has_access_to_genome($self->genome))) {
        CoGe::Exception::AccessDenied->throw(type => 'genome', id => $self->genome->id);
    }

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
