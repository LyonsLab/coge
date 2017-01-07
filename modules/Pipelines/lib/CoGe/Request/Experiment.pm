package CoGe::Request::Experiment;

use Moose;
#with qw(CoGe::Request::Request);
extends 'CoGe::Request::Request';

#use CoGe::Request::Request;
#use JSON;

sub is_valid {
    my $self = shift;

    # Verify that the experiment exists
    my $eid = $self->parameters->{eid} || $self->parameters->{experiment_id};
    return unless $eid;
    my $experiment = $self->db->resultset("Experiment")->find($eid);
    return defined $experiment ? 1 : 0;
}

sub has_access {
    my $self = shift;
    return unless defined $self->{user};

    my $eid = $self->parameters->{eid} || $self->parameters->{experiment_id};
    return unless $eid;
    my $experiment = $self->db->resultset("Experiment")->find($eid);
    return $self->user->has_access_to_experiment($experiment);
}

1;
