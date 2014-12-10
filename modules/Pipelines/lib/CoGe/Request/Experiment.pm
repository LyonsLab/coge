package CoGe::Requests::Experiment;

use Moose;
use CoGe::Requests::Request;
use JSON;

sub is_valid {
    my $self = shift;

    # Verify that the experiment exists
    my $eid = $self->parameters->{eid};
    my $experiment = $self->db->resultset("Experiment")->find($eid);
    return defined $experiment ? 1 : 0;
}

sub has_access {
    my $self = shift;

    my $eid = $self->parameters->{eid};
    my $experiment = $self->db->resultset("Experiment")->find($eid);
    return $self->user->has_access_to_experiment($experiment);
}

with qw(CoGe::Requests::Request);
1;
