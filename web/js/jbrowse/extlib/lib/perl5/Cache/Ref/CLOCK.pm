package Cache::Ref::CLOCK;
BEGIN {
  $Cache::Ref::CLOCK::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::CLOCK::VERSION = '0.04';
}
# ABSTRACT: CLOCK cache replacement algorithm

use Moose;

use namespace::autoclean;

extends qw(Cache::Ref);

with qw(Cache::Ref::CLOCK::Base);

has k => (
    isa => "Int",
    is  => "ro",
    default => 1,
);

sub _hit {
    my ( $self, $e ) = @_;

    my $k = 0+$self->k; # moose stingifies default
    $_->[0] = $k for @$e;
}

__PACKAGE__->meta->make_immutable;

__PACKAGE__;


# ex: set sw=4 et:

__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::CLOCK - CLOCK cache replacement algorithm

=head1 SYNOPSIS

    my $c = Cache::Ref::CLOCK->new(
        size => $n,
        k    => $k,
    );

=head1 DESCRIPTION

This algorithm is provides a second chance FIFO cache expiry policy using a
circular buffer.

It is a very well accepted page replacement algorithm, but largely for reasons
which are irrelevant in this context (cache hits don't need to be serialized in
a multiprocessing context as they only require an idempotent operation (setting
a bit to 1)).

=head1 ATTRIBUTES

=over 4

=item size

The size of the live entries.

=item k

This is the initial value given to all hit entries.

As the hand moves through the circular buffer it decrements the counters.

The default is C<1>, providing semantics similar to a second chance FIFO cache.

Larger values of C<k> model LRU more accurately.

This is pretty silly though, as L<Cache::Ref::LRU> is probably way more
efficient for any C<k> bigger than 1.

=back

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

