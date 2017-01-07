package CoGe::Request::Genome;

use Moose;
#with qw(CoGe::Request::Request);
extends 'CoGe::Request::Request';

sub is_valid {
    my $self = shift;

    my $gid = $self->parameters->{gid} || $self->parameters->{genome_id};
    unless ($gid) {
        warn "Request::Genome::is_valid: Missing gid/genome_id parameter";
        return;
    }

    my $genome = $self->db->resultset("Genome")->find($gid);
    unless ($genome) {
        warn "Request::Genome::is_valid: Genome $gid not found";
        return;
    }

    return 1;
}

sub has_access {
    my $self = shift;
    unless (defined $self->{user}) {
        warn "Request::Genome::has_access: User not defined";
        return;
    }

    my $gid = $self->parameters->{gid} || $self->parameters->{genome_id};
    unless ($gid) {
        warn "Request::Genome::has_access: Missing gid/genome_id parameter";
        return;
    }

    my $genome = $self->db->resultset("Genome")->find($gid);
    unless ($self->user->has_access_to_genome($genome)) {
        warn "Request::Genome::has_access: Access denied to genome $gid";
        return;
    }

    return 1;
}

1;
