package CoGe::Request::Empty;

use Moose;
#with qw(CoGe::Request::Request);
extends 'CoGe::Request::Request';

#use CoGe::Request::Request;
#use Data::Dumper;

sub is_valid {
#    my $self = shift;
    return 1;
}

sub has_access {
    my $self = shift;
    return unless defined $self->{user};

    return 1;
}

1;
