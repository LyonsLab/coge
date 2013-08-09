package CoGe::Accessory::Utils;

=head1 NAME

CoGe::Accessory::Utils

=head1 SYNOPSIS

Miscellaneous utility functions.

=head1 DESCRIPTION

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use POSIX qw( ceil );
use Sort::Versions;
use Data::Dumper;

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT =
      qw( units commify );
}

sub units {
    my $val = shift;

    if ( $val < 1024 ) {
        return $val;
    }
    elsif ( $val < 1024 * 1024 ) {
        return ceil( $val / 1024 ) . 'K';
    }
    elsif ( $val < 1024 * 1024 * 1024 ) {
        return ceil( $val / ( 1024 * 1024 ) ) . 'M';
    }
    else {
        return ceil( $val / ( 1024 * 1024 * 1024 ) ) . 'G';
    }
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

1;