package Cache::Ref::Util::LRU::API;
BEGIN {
  $Cache::Ref::Util::LRU::API::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Util::LRU::API::VERSION = '0.04';
}
use Moose::Role;

use namespace::autoclean;

requires qw(
    insert
    hit
    remove

    clear

    mru
    lru
    remove_mru
    remove_lru
);

__PACKAGE__;

# ex: set sw=4 et:


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Util::LRU::API

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

