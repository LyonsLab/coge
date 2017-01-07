package CoGe::Request::CoGeBlast;

use CoGe::Builder::Tools::CoGeBlast qw( get_genomes );
use Moose;
#with qw(CoGe::Request::Request);
extends 'CoGe::Request::Request';

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
    my @gids = get_genomes($self->parameters->{genomes}, $self->parameters->{notebooks}, $self->db);
    for (@gids) {
        if ($self->user) {
           return unless $self->user->has_access_to_genome($self->db->resultset("Genome")->find($_));
        }
        else {
            my @a = $self->db->storage->dbh->selectrow_array('SELECT restricted FROM genome WHERE genome_id=' . $_ );
            return if $a[0];
        }
    }
    return 1;
}

1;
