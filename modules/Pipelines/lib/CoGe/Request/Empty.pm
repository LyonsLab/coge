package CoGe::Request::Empty;

use Moose;
extends 'CoGe::Request::Request';

sub is_valid {
#    my $self = shift;
    return 1;
}

sub has_access {
    my $self = shift;
    return unless defined $self->{user};

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
