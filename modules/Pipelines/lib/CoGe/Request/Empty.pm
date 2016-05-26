package CoGe::Request::Empty;

use Moose;
with qw(CoGe::Request::Request);

use CoGe::Request::Request;

sub is_valid {
#    my $self = shift;
    return 1;
}

sub has_access {
    my $self = shift;
    return unless $self->parameters->{public} || defined $self->parameters->{user};

    return 1;
}

1;
