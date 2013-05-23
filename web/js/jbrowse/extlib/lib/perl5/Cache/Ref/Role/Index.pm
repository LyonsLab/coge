package Cache::Ref::Role::Index;
BEGIN {
  $Cache::Ref::Role::Index::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Role::Index::VERSION = '0.04';
}
use Moose::Role;

use namespace::autoclean;

# the index handles the by key lookup for all expiry methods
# the actual entries are set up by the manager though
# an entry in the index does not mean the key is live, it only means that it is
# known
has _index => (
    isa => "HashRef",
    default => sub { return {} },
    is => "ro",
);

sub _index_clear {
    my $self = shift;
    %{ $self->_index } = ();
}

sub _index_keys {
    my $self = shift;
    keys %{ $self->_index };
}

sub _index_get {
    my ( $self, @keys ) = @_;
    @{ $self->_index }{@keys};
}

sub _index_set {
    my ( $self, $key, $value ) = @_;
    $self->_index->{$key} = $value;
}

sub _index_size {
    my $self = shift;
    scalar keys %{ $self->_index };
}

sub _index_delete {
    my ( $self, @keys ) = @_;
    delete @{ $self->_index }{@keys};
}

# ex: set sw=4 et:

__PACKAGE__;


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Role::Index

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

