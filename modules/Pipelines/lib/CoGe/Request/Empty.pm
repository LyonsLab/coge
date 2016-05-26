package CoGe::Request::Empty;

use Moose;
with qw(CoGe::Request::Request);

use CoGe::Request::Request;
use Data::Dumper;

sub is_valid {
#    my $self = shift;
    return 1;
}

sub has_access {
    my $self = shift;
    warn Dumper $self->options;
    warn Dumper $self->parameters;
    return unless ($self->options && $self->options->{public}) || defined $self->parameters->{user};

    return 1;
}

1;
