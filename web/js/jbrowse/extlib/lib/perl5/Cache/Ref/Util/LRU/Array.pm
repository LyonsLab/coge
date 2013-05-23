package Cache::Ref::Util::LRU::Array;
BEGIN {
  $Cache::Ref::Util::LRU::Array::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Util::LRU::Array::VERSION = '0.04';
}
use Moose;

use Scalar::Util qw(refaddr);
use Hash::Util::FieldHash::Compat qw(id);

use namespace::autoclean;

has _list => (
    traits => [qw(Array)],
    isa => "ArrayRef",
    default => sub { [] },
    is => "ro",
    handles => {
        #size       => "length",
        mru        => [ get => 0 ],
        lru        => [ get => -1 ],
        remove_mru => "shift",
        remove_lru => "pop",
        clear      => "clear",
    },
);

with qw(Cache::Ref::Util::LRU::API);

# since there's no need for metadata, insert is just like hit
sub insert {
    my ( $self, @elements ) = @_;

    $self->hit(@elements);

    return ( @elements == 1 ? $elements[0] : @elements );
}

sub _filter {
    my ( $self, $l, $elements ) = @_;

    return () unless @$l;

    confess if grep { not defined } @$elements;
    my %hash; @hash{map {id($_)} @$elements} = ();
    grep { not exists $hash{id($_)} } @$l;
}

sub hit {
    my ( $self, @elements ) = @_;

    return unless @elements;

    my $l = $self->_list;
    @$l = ( @elements, $self->_filter($l, \@elements) );

    return;
}

sub remove {
    my ( $self, @elements ) = @_;

    return unless @elements;

    my $l = $self->_list;
    @$l = $self->_filter($l, \@elements);

    return;
}

__PACKAGE__->meta->make_immutable;

__PACKAGE__;

# ex: set sw=4 et:


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Util::LRU::Array

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

