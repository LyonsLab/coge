package Cache::Ref::LRU;
BEGIN {
  $Cache::Ref::LRU::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::LRU::VERSION = '0.04';
}
# ABSTRACT: Least recently used expiry policy

use Moose;

use Cache::Ref::Util::LRU::List;

use namespace::autoclean;

extends qw(Cache::Ref);

with qw(
    Cache::Ref::Role::API
    Cache::Ref::Role::Index
);

has size => (
    isa => "Int",
    is  => "ro",
    required => 1,
);

has lru_class => (
    isa => "ClassName",
    is  => "ro",
    default => "Cache::Ref::Util::LRU::List",
);

has _lru => (
    does => "Cache::Ref::Util::LRU::API",
    is   => "ro",
    lazy_build => 1,
);

sub _build__lru { shift->lru_class->new }

sub get {
    my ( $self, @keys ) = @_;

    my @e = $self->_index_get(@keys);

    $self->_lru->hit(map { $_->[1] } grep { defined } @e);

    return ( @keys == 1 ? $e[0][0] : map { $_ && $_->[0] } @e );
}

sub hit {
    my ( $self, @keys ) = @_;

    $self->_lru->hit( map { $_->[1] } $self->_index_get(@keys) );

    return;
}

sub expire {
    my ( $self, $how_many ) = @_;

    my $l = $self->_lru;
    $self->_index_delete( $l->remove_lru ) for 1 .. ($how_many || 1);

    return;
}

sub set {
    my ( $self, $key, $value ) = @_;

    my $l = $self->_lru;

    if ( my $e = $self->_index_get($key) ) {
        $l->hit($e->[1]);
        $e->[0] = $value;
    } else {
        if ( $self->_index_size == $self->size ) {
            $self->expire(1);
        }

        $self->_index_set( $key => [ $value, $l->insert($key) ] );
    }

    return $value;
}

sub clear {
    my $self = shift;

    $self->_lru->clear;
    $self->_index_clear;

    return;
}

sub remove {
    my ( $self, @keys ) = @_;

    $self->_lru->remove(map { $_->[1] } $self->_index_delete(@keys));

    return;
}

__PACKAGE__->meta->make_immutable;

__PACKAGE__;


# ex: set sw=4 et:

__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::LRU - Least recently used expiry policy

=head1 SYNOPSIS

    my $c = Cache::Ref::LRU->new(
        size => $n,
    );

=head1 DESCRIPTION

This is an implementation of the least recently used expiry policy.

It provides both an array and a doubly linked list based implementation. See
L<Cache::Ref> for a discussion.

=head1 ATTRIBUTES

=over 4

=item size

The size of the live entries.

=item lru_class

The class of the LRU list implementation.

=back

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

