package Data::Stag::Arr2HTML;

=head1 NAME

  Data::Stag::Arr2HTML

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Base);

use vars qw($VERSION);
$VERSION="0.11";

sub start_event {
    my $self = shift;
    my $ev = shift;
    print "<tr><td>$ev<td>\n";
}

sub end_event {
    my $self = shift;
    my $ev = shift;
    print "</tr>\n";
}

sub evbody {
    my $self = shift;
    my $body = shift;
    print "<td>$body</td>\n";
}

1;
