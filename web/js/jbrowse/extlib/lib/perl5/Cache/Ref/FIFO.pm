package Cache::Ref::FIFO;
BEGIN {
  $Cache::Ref::FIFO::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::FIFO::VERSION = '0.04';
}
use Moose;

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

has _fifo => (
    isa => "ArrayRef",
    is  => "ro",
    lazy => 1,
    default => sub { [] },
);

sub clear {
    my $self = shift;
    $self->_index_clear;
    @{ $self->_fifo } = ();
}

sub hit { }

sub remove {
    my ( $self, @keys ) = @_;

    $self->_index_delete(@keys);

    my %keys; @keys{@keys} = ();

    my $f = $self->_fifo;
    @$f = grep { not exists $keys{$_} } @$f;

    return;
}

sub get {
    my ( $self, @keys ) = @_;
    $self->_index_get(@keys);
}

sub set {
    my ( $self, $key, $value ) = @_;

    unless ( defined $self->_index_get($key) ) {
        if ( $self->_index_size >= $self->size ) {
            $self->expire( 1 + $self->_index_size - $self->size );
        }
        push @{ $self->_fifo }, $key;
    }

    $self->_index_set($key, $value);
}

sub expire {
    my ( $self, $how_many ) = @_;

    $self->_index_delete( splice @{ $self->_fifo }, 0, $how_many || 1 );

    return;
}

__PACKAGE__->meta->make_immutable;

# ex: set sw=4 et:

__PACKAGE__;


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::FIFO

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

