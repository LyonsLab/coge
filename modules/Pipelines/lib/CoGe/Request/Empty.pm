package CoGe::Request::Empty;

use Moose;
use CoGe::Request::Request;

sub is_valid {
    my $self = shift;
    return 1;
}

sub has_access {
    my $self = shift;
    return 1;
}

with qw(CoGe::Request::Request);
1;
