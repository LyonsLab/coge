package CoGe::Request::ExperimentAnalysis;

use Moose;
extends 'CoGe::Request::Request';

has genome     => (is => 'rw', isa => 'CoGeX::Result::Genome');
has experiment => (is => 'rw', isa => 'CoGeX::Result::Experiment');

sub is_valid { # called first
    my $self = shift;

    # Verify that the experiment exists
    my $eid = $self->parameters->{eid} || $self->parameters->{experiment_id};
    return unless $eid;

    my $experiment = $self->db->resultset("Experiment")->find($eid);
    if ($experiment) {
        $self->experiment( $experiment );
        $self->genome( $experiment->genome );
        return 1;
    }

    return;
}

sub has_access { # called second
    my $self = shift;
    return unless defined $self->{user};

    my $eid = $self->parameters->{eid} || $self->parameters->{experiment_id};
    return unless $eid;
    my $experiment = $self->db->resultset("Experiment")->find($eid);
    return $self->user->has_access_to_genome($experiment->genome);
}

__PACKAGE__->meta->make_immutable;

1;
