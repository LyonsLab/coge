package CoGe::Request::NExperiments;

use Moose;
extends 'CoGe::Request::Request';

use CoGe::Exception::Generic;
use CoGe::Exception::AccessDenied;
use CoGe::Exception::ItemNotFound;

has experiments => (is => 'rw', isa => 'ArrayRef', default => sub{ [] });

sub is_valid {
    my $self = shift;

    unless ($self->parameters->{experiments}) {
        CoGe::Exception::Generic->throw(message => 'Missing required "experiments" parameter');
    }

	foreach my $eid (@{$self->parameters->{experiments}}) {
        my $experiment = $self->db->resultset("Experiment")->find($eid);
        unless ($experiment) {
            CoGe::Exception::ItemNotFound->throw(type => 'experiment', id => $eid);
        }

        push @{$self->experiments}, $experiment;
	}

    if (scalar(@{$self->experiments}) < 2) {
        CoGe::Exception::Generic->throw(message => 'At least 2 experiments must be specified');
    }

    return 1;
}

sub has_access {
    my $self = shift;

    foreach my $experiment (@{$self->experiments}) {
        if ($experiment->restricted && (!$self->user || !$self->user->has_access_to_experiment($experiment))) {
            CoGe::Exception::AccessDenied->throw(type => 'experiment', id => $experiment->id);
        }
    }

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
