package Cache::Ref::Util::LRU::List;
BEGIN {
  $Cache::Ref::Util::LRU::List::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Util::LRU::List::VERSION = '0.04';
}
use Moose;

use namespace::autoclean;

with (
    'Cache::Ref::Util::LRU::API',
    'Cache::Ref::Role::WithDoublyLinkedList' => {
        name => "",
        value_offset => 0,
        next_offset => 1,
        prev_offset => 2,

        head_method  => "mru",
        tail_method  => "lru",
        shift_method => "remove_mru",
        pop_method   => "remove_lru",
        unshift_method => "insert",
    },
);

sub hit {
    my ( $self, @l ) = @_;

    return unless @l;

    $self->_splice(@l);
    $self->_link_sequence(@l);

    my ( $head, $tail ) = @l[0, -1];

    if ( my $old_head = $self->_head ) {
        # prepend @l to the head of the list
        $self->_set_prev($old_head, $tail);
        $self->_set_next($tail, $old_head);
    } else {
        # just set @l as the new list
        $self->_tail($tail);
    }

    $self->_head($head);

    return;
}

sub remove {
    my ( $self, @l ) = @_;

    $self->_splice(@l);
    @$_ = () for @l;
}

sub clear { shift->_clear }

sub DEMOLISH { shift->clear }

__PACKAGE__->meta->make_immutable;

__PACKAGE__;

# ex: set sw=4 et:


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Util::LRU::List

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

