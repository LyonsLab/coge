package CoGe::Request::CoGeBlast;

use CoGe::Builder::Tools::CoGeBlast qw( get_genomes );
use Moose;
with qw(CoGe::Request::Request);

use CoGe::Request::Request;
use JSON;

sub is_valid {
    my $self = shift;
    return unless $self->parameters->{query_seq};
    return unless scalar get_genomes($self->parameters->{genomes}, $self->parameters->{notebooks}, $self->db);
    return 1;
}

sub has_access {
    my $self = shift;
    return unless defined $self->{user};
    my @gids = get_genomes($self->parameters->{genomes}, $self->parameters->{notebooks}, $self->db);
    for (@gids) {
        return unless $self->user->has_access_to_genome($self->db->resultset("Genome")->find($_));
    }
    return 1;
}

1;
