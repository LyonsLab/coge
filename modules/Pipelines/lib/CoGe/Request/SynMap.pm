package CoGe::Request::SynMap;

use Moose;
with qw(CoGe::Request::Request);

use CoGe::Request::Request;
use Data::Dumper;

sub is_valid {
#    my $self = shift;
    return 1;
}

sub has_access {
#    my $self = shift;
    return 1;
}

1;
