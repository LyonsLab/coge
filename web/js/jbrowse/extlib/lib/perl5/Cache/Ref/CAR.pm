package Cache::Ref::CAR;
BEGIN {
  $Cache::Ref::CAR::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::CAR::VERSION = '0.04';
}
# ABSTRACT: CLOCK with Adaptive Replacement

use Moose;

use List::Util qw(max min);
use Cache::Ref::CAR::Base ();

use namespace::autoclean;

extends qw(Cache::Ref);

with qw(Cache::Ref::CAR::Base);

sub _mru_history_too_big {
    my $self = shift;

    $self->_mru_history_size
        and
    $self->_mru_history_size + $self->_mru_size == $self->size;
}

sub _mfu_history_too_big {
    my $self = shift;

    $self->_index_size == $self->size * 2;
}

sub _decrease_mru_target_size {
    my $self = shift;

    if ( $self->_mru_target_size > 0 ) {
        my $adjustment = int( $self->_mru_history_size / $self->_mfu_history_size );
        $self->_set_mru_target_size( max( 0, $self->_mru_target_size - max(1, $adjustment) ) );
    }
}

sub _increase_mru_target_size {
    my $self = shift;

    my $adjustment = int( $self->_mfu_history_size / $self->_mru_history_size );
    $self->_set_mru_target_size( min( $self->size, $self->_mru_target_size + max(1, $adjustment) ) );
}

sub _restore_from_mfu_history {
    my ( $self, $e ) = @_;

    $self->_mfu_push($e);
}

sub _restore_from_mru_history {
    my ( $self, $e ) = @_;

    $self->_mfu_push($e);
}

sub expire {
    my ( $self, $how_many ) = @_;

    $how_many ||= 1;

    if ( my $mru = $self->_mru ) {
        my $cur = $self->_next($mru);

        # mru pool is too big
        while ( $cur and $self->_mru_size >= max(1,$self->_mru_target_size) ) {
            my $next = $self->_next($cur);

            $self->_circular_splice($cur);

            if ( $cur->[0] & Cache::Ref::CAR::Base::REF_BIT ) {
                $cur->[0] &= ~Cache::Ref::CAR::Base::REF_BIT; # turn off reference bit

                # move to t2 (mfu)
                $self->_mfu_push($cur);

                $cur = $next;
            } else {
                # reference bit is off, which means this entry is freeable

                delete $cur->[2]; # delete the value

                # move to history
                # MFU_BIT not set

                if ( $self->_mru_history_head ) {
                    $self->_set_next($cur, $self->_mru_history_head);
                    $self->_set_prev($self->_mru_history_head, $cur);
                } else {
                    $self->_set_next($cur, undef);
                }

                $self->_mru_history_head($cur);
                $self->_mru_history_tail($cur) unless $self->_mru_history_tail;
                $self->_inc_mru_history_size;


                return;
            }
        }
    }

    for ( 1 .. $how_many ) {
        my $tail = $self->_mfu || last;
        my $cur = $self->_next($tail) || last;

        loop: {
            if ( $cur->[0] & Cache::Ref::CAR::Base::REF_BIT ) {
                $cur->[0] &= ~Cache::Ref::CAR::Base::REF_BIT;
                $tail = $cur;
                $cur = $self->_next($cur);
                redo loop;
            } else {
                # reference bit is off, which means this entry is freeable

                $self->_mfu($tail);
                $self->_circular_splice($cur);

                delete $cur->[2]; # delete the value

                # move to history
                $cur->[0] |= Cache::Ref::CAR::Base::MFU_BIT;

                if ( $self->_mfu_history_head ) {
                    $self->_set_prev($self->_mfu_history_head, $cur);
                    $self->_set_next($cur, $self->_mfu_history_head);
                } else {
                    $self->_set_next($cur, undef);
                }

                $self->_mfu_history_head($cur);
                $self->_mfu_history_tail($cur) unless $self->_mfu_history_tail;
                $self->_inc_mfu_history_size;
            }
        }
    }

    return;
}

sub _clear_additional { }

__PACKAGE__->meta->make_immutable;

__PACKAGE__;


# ex: set sw=4 et:


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::CAR - CLOCK with Adaptive Replacement

=head1 SYNOPSIS

    my $c = Cache::Ref::CAR->new(
        size => $n,
    );

=head1 DESCRIPTION

This algorithm is an implementation of
L<http://www.almaden.ibm.com/cs/people/dmodha/clockfast.pdf|CAR: Clock with Adaptive Replacement, Sorav Bansal and Dharmendra S. Modha>.

See also L<Cache::Ref::CART> which is probably more appropriate for random access work loads.

CAR balances between an MFU like policy and an MRU like policy, automatically
tuning itself as the workload varies.

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

