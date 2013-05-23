package Cache::Ref::CART;
BEGIN {
  $Cache::Ref::CART::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::CART::VERSION = '0.04';
}
# ABSTRACT: CAR with temporal filtering

use Moose;

use List::Util qw(max min);
use Cache::Ref::CAR::Base ();

use namespace::autoclean;

extends qw(Cache::Ref);

with qw(Cache::Ref::CAR::Base);

has _long_term_utility_size => ( # q
    is => "ro",
    writer => "_set_long_term_utility_size",
    default => 0,
);

has _mru_history_target_size => ( # q
    is => "ro",
    writer => "_set_mru_history_target_size",
    default => 0,
);

sub _inc_long_term_utility_size {
    my $self = shift;
    $self->_set_long_term_utility_size( $self->_long_term_utility_size + 1 );
}

sub _dec_long_term_utility_size {
    my $self = shift;
    $self->_set_long_term_utility_size( $self->_long_term_utility_size - 1 );
}

sub _reset_long_term_utility_size {
    my $self = shift;
    $self->_set_long_term_utility_size(0);
}

sub _reset_mru_history_target_size {
    my $self = shift;
    $self->_set_mru_history_target_size(0);
}

sub _mru_history_too_big {
    my $self = shift;

    # only if there's something to purge
    return unless $self->_mru_history_size;

    # only if we need to purge
    return unless $self->_mru_history_size + $self->_mfu_history_size == $self->size + 1;

    # purge from here if there's nothing to purge from mfu
    return 1 if $self->_mfu_history_size == 0;

    # or if the target size is too big
    return 1 if $self->_mru_history_size > $self->_mru_history_target_size;

    return;
}

sub _mfu_history_too_big {
    my $self = shift;

    return unless $self->_mfu_history_size;

    # only purge if we actually need to
    return 1 if $self->_mru_history_size + $self->_mfu_history_size == $self->size + 1;

    return;
}

sub _increase_mru_target_size {
    my $self = shift;

    my $adjustment = int( ( $self->_mru_size + $self->_mfu_size - $self->_long_term_utility_size ) / $self->_mru_history_size );
    $self->_set_mru_target_size( min( $self->size, $self->_mru_target_size + max(1, $adjustment) ) );
}

sub _decrease_mru_target_size {
    my $self = shift;

    my $adjustment = int( $self->_long_term_utility_size / $self->_mfu_history_size );
    $self->_set_mru_target_size( max( 0, $self->_mru_target_size - max(1, $adjustment) ) );
}

sub _increase_mru_history_target_size {
    my $self = shift;

    $self->_set_mru_history_target_size( min($self->_mru_history_target_size + 1, 2 * $self->size - $self->_mru_size ) );
}

sub _decrease_mru_history_target_size {
    my $self = shift;

    $self->_set_mru_history_target_size( max($self->_mru_history_target_size - 1, $self->size - $self->_mru_size) );
}

sub _restore_from_mfu_history {
    my ( $self, $e ) = @_;

    if ( $self->_mfu_history_size + $self->_long_term_utility_size >= $self->size ) {
        $self->_increase_mru_history_target_size();
    }

    die unless $e->[0] & Cache::Ref::CAR::Base::LONG_TERM_BIT();
    $self->_inc_long_term_utility_size();

    $self->_mru_push($e);
}

sub _restore_from_mru_history {
    my ( $self, $e ) = @_;

    $e->[0] |= Cache::Ref::CAR::Base::LONG_TERM_BIT();
    $self->_inc_long_term_utility_size();

    $self->_mru_push($e);
}

sub expire {
    my ( $self, $how_many ) = @_;

    $how_many ||= 1;

    if ( my $mfu = $self->_mfu ) {
        my $cur = $self->_next($mfu);

        # mru pool is too big
        while ( $cur and $cur->[0] & Cache::Ref::CAR::Base::REF_BIT ) {
            $self->_circular_splice($cur);

            $cur->[0] &= ~Cache::Ref::CAR::Base::REF_BIT; # turn off reference bit

            # move to t1 (mru)
            $self->_mru_push($cur);
            $cur = $self->_next($cur);

            if ( $self->_mfu_history_size + $self->_long_term_utility_size >= $self->size ) {
                $self->_increase_mru_history_target_size;
            }
        }
    }

    if ( my $mru = $self->_mru ) {
        my $cur = $self->_next($mru);

        while ( $cur ) {
            if ( $cur->[0] & Cache::Ref::CAR::Base::REF_BIT ) {
                $cur->[0] &= ~Cache::Ref::CAR::Base::REF_BIT;

                if ( $self->_mru_size >= max($self->_mru_history_size, $self->_mru_target_size + 1)
                        and
                    not( $cur->[0] & Cache::Ref::CAR::Base::LONG_TERM_BIT )
                ) {
                    # FIXME spec says 'x', is this the same as 'head'?
                    $cur->[0] |= Cache::Ref::CAR::Base::LONG_TERM_BIT;
                    $self->_inc_long_term_utility_size;
                }

                # $hand++
                $self->_mru($cur);
                $cur = $self->_next($cur);
            } elsif ( $cur->[0] & Cache::Ref::CAR::Base::LONG_TERM_BIT ) {
                my $next = $self->_next($cur);
                $self->_circular_splice($cur);
                $self->_mfu_push($cur);
                $cur = $self->_next($self->_mru);;

                $self->_decrease_mru_history_target_size();
            } else {
                # found a candidate page for removal
                last;
            }
        }
    }

    for ( 1 .. $how_many ) {
        if ( $self->_mru_size >= max(1, $self->_mru_target_size) ) {
            my $head = $self->_next($self->_mru);
            $self->_circular_splice($head);

            if ( $self->_mru_history_head ) {
                $self->_set_next($head, $self->_mru_history_head);
                $self->_set_prev($self->_mru_history_head, $head);
            }

            $self->_mru_history_head($head);
            $self->_mru_history_tail($head) unless $self->_mru_history_tail;
            $self->_inc_mru_history_size;

            delete $head->[2]; # delete the value
        } else {
            my $tail = $self->_mfu || last;
            my $head = $self->_next($tail) || last;
            $self->_circular_splice($head);

            $self->_dec_long_term_utility_size; # entries in mfu *always* have long term set

            if ( $self->_mfu_history_head ) {
                $self->_set_next($head, $self->_mfu_history_head);
                $self->_set_prev($self->_mfu_history_head, $head);
            }

            $self->_mfu_history_head($head);
            $self->_mfu_history_tail($head) unless $self->_mfu_history_tail;
            $self->_inc_mfu_history_size;

            delete $head->[2]; # delete the value
            $head->[0] |= Cache::Ref::CAR::Base::MFU_BIT;
        }
    }

    return;
}

sub _clear_additional {
    my $self = shift;

    $self->_reset_long_term_utility_size;
    $self->_reset_mru_history_target_size;
}

__PACKAGE__->meta->make_immutable;

__PACKAGE__;


# ex: set sw=4 et:


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::CART - CAR with temporal filtering

=head1 SYNOPSIS

    my $c = Cache::Ref::CART->new(
        size => $n,
    );

=head1 DESCRIPTION

This algorithm is an extension to L<Cache::Ref::CAR> that has temporal
filtering on the upgrading from MRU to MFU pool.

This means that two subsequent accesses to the same key do not automatically
make it viable for long term caching, to get upgraded to MFU status a key must
be expired but known in the history.

This is probably the most general purpose caching algorithm.

=head1 ATTRIBUTES

=over 4

=item size

The size of the live entries.

Note that the cache also remembers this many expired keys, and keeps some
metadata about those keys, so for memory usage the overhead is probably around
double what L<Cache::Ref::LRU> requires.

=back

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

