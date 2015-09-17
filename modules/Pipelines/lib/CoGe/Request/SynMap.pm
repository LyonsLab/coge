package CoGe::Request::Empty;

use Moose;
use CoGe::Request::Request;

sub is_valid {
    my $self = shift;
    my $gid1 = $self->parameters->{gid1};
    return unless $gid1;
    my $gid2 = $self->parameters->{gid2};
    return unless $gid2;
    my $genome1 = $self->db->resultset("Genome")->find($gid1);
    return unless defined $genome1;
    my $genome2 = $self->db->resultset("Genome")->find($gid2);
    return unless defined $genome2;
    return 1;
}

sub has_access {
    my $self = shift;
    my $gid1 = $self->parameters->{gid1};
    return unless $gid1;
    my $gid2 = $self->parameters->{gid2};
    return unless $gid2;
    my $genome = $self->db->resultset("Genome")->find($gid);
    return $self->user->has_access_to_genome($genome);
    return 1;
}

with qw(CoGe::Request::Request);
1;
