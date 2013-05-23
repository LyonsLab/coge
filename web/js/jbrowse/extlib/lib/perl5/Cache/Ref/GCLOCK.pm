package Cache::Ref::GCLOCK;
BEGIN {
  $Cache::Ref::GCLOCK::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::GCLOCK::VERSION = '0.04';
}
# ABSTRACT: GCLOCK cache replacement algorithm

use Moose;

use namespace::autoclean;

extends qw(Cache::Ref);

with qw(Cache::Ref::CLOCK::Base);

sub _hit {
    my ( $self, $e ) = @_;

    $_->[0]++ for @$e;
}

__PACKAGE__->meta->make_immutable;

__PACKAGE__;


# ex: set sw=4 et:


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::GCLOCK - GCLOCK cache replacement algorithm

=head1 SYNOPSIS

    my $c = Cache::Ref::GCLOCK->new(
        size => $n,
    );

=head1 DESCRIPTION

This algorithm is related to L<Cache::Ref::CLOCK> but instead of starting all
cache hits from C<k>, a counter is increased on every hit.

This provides behavior which models an LFU expiry policy (without taking into
account the full keyspace).

=head1 ATTRIBUTES

=over 4

=item size

The size of the live entries.

=back

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

