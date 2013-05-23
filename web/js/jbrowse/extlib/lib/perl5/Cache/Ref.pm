package Cache::Ref;
BEGIN {
  $Cache::Ref::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::VERSION = '0.04';
}
# ABSTRACT: Memory only cache of live references

use Moose;

__PACKAGE__->meta->make_immutable;

__PACKAGE__;


# ex: set sw=4 et:

__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref - Memory only cache of live references

=head1 SYNOPSIS

    # this class is just a base class and a documentation start point
    # just use the various algorithms directly

    use Cache::Ref::CART;
    my $cache = Cache::Ref::CART->new( size => 1024 );


    # add a cache value or set an existing key to a new value
    $cache->set(foo => $some_object);


    # get a value
    $cache->get("foo"); # also takes a list of keys


    # remove a key before it has normally expired
    $cache->remove("foo");


    # remove all cached data
    $cache->clear;


    # 'hit' is like 'get' without the overhead of obtaining the value
    # it's useful for keeping values from expiring when you already have
    # the values
    $cache->hit("foo"); # also takes a list of keys

=head1 DESCRIPTION

Unlike L<CHI> which attempts to address the problem of caching things
persistently, this module implements in memory caching, designed primarily for
B<shared references> in memory.

This collection of classes implements a number of semi related algorithms.

=head1 METHODS

=over 4

=item get @keys

Fetch entries from the cache.

=item hit @keys

Promote C<@keys> in the cache.

Same effect as C<get> except it doesn't actually return anything.

=item set $key => $value

Adds an entry to the cache.

=item compute $key, sub { ...; return $value }

Calls C<get> with C<$key>. If there's a hit the value is
returned. Otherwise the code block is executed to compute the value, and the result is stored in the cache using C<set>.

=item remove @keys

Remove specific entries from the cache.

=item expire $x

Remove C<$x> many entries from the cache. Hopefully the entries
removed are the most useless ones.

C<$x> defaults to 1.

=item clear

Empty the cache.

=back

=head1 ALGORITHMS

=head2 FIFO

This is a simple FIFO queue where a C<set> places the element on the head of a
queue, and if the size is too big an element will be discarded from the tail of
the queue.

L<Cache::Bounded> provides similar behavior, but flushing happens periodically
and in bigger numbers. Therefore, performance will be better on very high cache
usage, when hits don't matter that much.

This implementation has the lowest memory overhead, due to the simplicity of
its data structures (just a hash and an array).

Its expiry policy is appropriate for when the data set has a high locality of
reference, and random access is generally confined to neighbors, as a part of
some larger scan.

For truly random access cache hit rates will suffer.

Long term utility of cache entries is not considered at all, so scans will
poison the cache.

This is the only algorithm for which C<get> (and C<hit>) has no side effects.

=head2 LRU

This implementation uses an LRU list of entries (two implementations are
provided for trading off memory for speed).

Long term utility of cache entries is not considered at all, so scans will
poison the cache.

=head3 Cache::Ref::Util::LRU::List

Uses a doubly linked list to perform MRU propagation.

Faster than Array.

Cache hits and LRU removal is O(1).

=head3 Cache::Ref::Util::LRU::Array

Generally slower for a cache size bigger than about 10 elements, but uses less memory due to the compact layout.

Cache hits are O(cache size). LRU removal is O(1).

=head2 CLOCK

This is an implementation of second chance FIFO, using a circular buffer.

Second chance FIFO is a very simple approximation of LRU. The CLOCK algorithm
has its origins in Multics' virtual memory paging implementation.

It's slightly more general purpose than FIFO when dealing with random access.

Long term utility of cache entries is not considered at all, so scans will
poison the cache.

Using values of C<k> bigger than 1 (the default), more accurate approximations
of LRU can be made, at the cost of more complicated expiry.

=head2 GCLOCK

Tries to approximate LFU instead of LRU.

Cache hits increment a counter by one, instead of resetting it to the constant C<k>.

Cache replacement decays existing counters just like CLOCK.

=head2 CAR

CLOCK with Adaptive Removal.

A self tuning cache that varies between approximations of LRU and LFU expiry.

Has the highest memory overhead of all the implementations due to the extent of
the metadata it maintains.

However, this overhead is still small for when sizeable objects are involved.

Resistent to cache poisoning when scanning.

=head2 CART

CAR with temporal filtering.

Like CAR but does not promote a cache entry to the long term usefulness set due
to frequent successive access.

This is probably the most general purpose algorithm.

=head1 SEE ALSO

=over 4

=item L<CHI>

Appropriate for persistent caching of data with complex expiry.

=item L<Cache::Cascade>

Can be used to layer L<Cache::Ref> over other caches (e.g. L<CHI>).

=item L<Cache::Bounded>

A simpler implementation with similar goals (memory only caching), designed for
when cache misses are not very high cost, so cache hits have an extremely low
overhead and the policy is very simplistic.

=item L<Cache::Weak>

Caches shared references for as long as there is some other reference to those
objects.

=item L<Cache::Profile>

Designed to help choose an appropriate cache layer.

=item Algorithm information

L<http://en.wikipedia.org/wiki/Cache_algorithms>

L<http://en.wikipedia.org/wiki/Page_replacement_algorithm>

L<http://www.almaden.ibm.com/cs/people/dmodha/clockfast.pdf>

=back

=head1 VERSION CONTROL

L<http://github.com/nothingmuch/Cache-Ref>

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

