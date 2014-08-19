package CoGe::Core::List;

use v5.14;
use strict;
use warnings;

our (@ISA, @EXPORT_OK, $VERSION);

BEGIN {
    require Exporter;
    $VERSION = 0.1;

    @ISA = qw(Exporter);
    @EXPORT_OK = qw(listcmp);
}

sub listcmp($$) {
    my ($a, $b) = $_;

    my $namea = "";
    my $nameb = "";

    $namea = $a->name if defined $a and $a->name;
    $nameb = $b->name if defined $b and $b->name;

    $namea cmp $nameb;
}

1;
