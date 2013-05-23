package Cache::Ref::CLOCK::Base;
BEGIN {
  $Cache::Ref::CLOCK::Base::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::CLOCK::Base::VERSION = '0.04';
}
use Moose::Role;

use namespace::autoclean;

with qw(
    Cache::Ref::Role::API
    Cache::Ref::Role::Index
);

requires qw(_hit);

has size => (
    isa => "Int",
    is  => "ro",
    required => 1,
);

has _hand => (
    isa => "ScalarRef",
    is  => "ro",
    default => sub { my $x = -1; return \$x },
);

has _buffer => (
    isa => "ArrayRef",
    is  => "ro",
    lazy => 1,
    default => sub {
        my ( $self, $p ) = @_;
        return [ map { [] } 1 .. $self->size ],
    },
);

sub hit {
    my ( $self, @keys ) = @_;

    $self->_hit( [ grep { defined } $self->_index_get(@keys) ] );

    return;
}

sub get {
    my ( $self, @keys ) = @_;

    my @ret;

    my @entries = $self->_index_get(@keys);

    $self->_hit( [ grep { defined } @entries ] );

    return ( @keys == 1 ? ($entries[0] && $entries[0][2]) : map { $_ && $_->[2] } @entries );
}

sub clear {
    my $self = shift;
    $self->_index_clear;
    @$_ = () for @{ $self->_buffer };
}

sub remove {
    my ( $self, @keys ) = @_;
    @$_ = () for $self->_index_delete(@keys);
}

sub set {
    my ( $self, $key, $value ) = @_;

    if ( my $e = $self->_index_get($key) ) {
        $e->[2] = $value;
    } else {
        my $e = $self->_find_free_slot;
        @$e = ( 0, $key, $value ); # start at 0, not k
        $self->_index_set( $key, $e );
    }
}

sub expire {
    my ( $self, $how_many ) = @_;

    my $i = $self->_hand;
    my $b = $self->_buffer;

    while ( $how_many ) {
        if ( $$i == $#$b ) {
            $$i = -1;
        }

        if ( my $e = $b->[++$$i] ) {
            if ( !$e->[0] ) {
                $self->remove($e->[1]); # also clears @$e
                $how_many--;
            } else {
                $e->[0]--;
            }
        }
    }

    return;
}

sub _find_free_slot {
    my $self = shift;

    my $i = $self->_hand;
    my $b = $self->_buffer;

    loop: {
        if ( $$i == $#$b ) {
            $$i = -1;
        }

        my $e = $b->[++$$i];

        if ( not @$e ) {
            return $e;
        } elsif ( !$e->[0] ) {
            $self->remove($e->[1]); # also clears @$e
            return $e;
        } else {
            $e->[0]--;
            redo loop;
        }
    }
}

# ex: set sw=4 et:

__PACKAGE__;


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::CLOCK::Base

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

