package Cache::Ref::LIFO;
BEGIN {
  $Cache::Ref::LIFO::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::LIFO::VERSION = '0.04';
}
# ABSTRACT: Saves entries until full, discarding subsequent sets.

use Moose;

use namespace::autoclean;

extends qw(Cache::Ref);

with (
    'Cache::Ref::Role::API',
    'Cache::Ref::Role::Index' => {
        -alias => {
            _index_get => "get",
        },
    },
);

has size => (
    isa => "Int",
    is  => "ro",
    required => 1,
);

has _keys => (
    traits => [qw(Array)],
    isa => "ArrayRef",
    is  => "ro",
    default => sub { [] },
    handles => {
        _add_key => "push",
        #_splice_keys => "splice",
        _clear_keys => "clear",
    },
);

sub remove {
    my ( $self, @keys ) = @_;

    $self->_index_delete(@keys);

    my %keys;
    undef @keys{@keys};

    @{ $self->_keys } = grep { not exists $keys{$_} } @{ $self->_keys };

    return;
}

sub clear {
    my $self = shift;
    $self->_index_clear;
    $self->_clear_keys;

    return;
}

sub hit { }

sub set {
    my ( $self, $key, $value ) = @_;

    if ( defined $self->_index_get($key) ) {
        $self->_index_set($key, $value)
    } elsif ( $self->_index_size < $self->size ) {
        $self->_index_set($key, $value);
        $self->_add_key($key);
    }

    return $value;
}

sub expire {
    my ( $self, $how_many ) = @_;

    $how_many ||= 1;

    #my @keys = $self->_splice_keys( -$how_many );
    my @keys = splice @{ $self->_keys }, -$how_many;

    $self->_index_delete(@keys);
}

__PACKAGE__->meta->make_immutable;

__PACKAGE__;



=pod

=encoding utf-8

=head1 NAME

Cache::Ref::LIFO - Saves entries until full, discarding subsequent sets.

=head1 SYNOPSIS

    my $c = Cache::Ref::LIFO->new( size => $n );

    $c->set( foo => 42 );

    $c->get("foo");

=head1 DESCRIPTION

This is a very naive cache algorithm, it saves cache sets until the
cache is full, at which point all additional saves which aren't a
value update are discarded immediately.

For very predictable workflows this is potentially a good fit,
provided the MFU is used early on.

The advantages is that the code is very simple as a result.

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut


__END__


# ex: set sw=4 et:
