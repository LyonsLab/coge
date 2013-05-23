package Cache::Ref::Null;
BEGIN {
  $Cache::Ref::Null::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Null::VERSION = '0.04';
}
# ABSTRACT: Caches nothing

use Moose;

use namespace::autoclean;

extends qw(Cache::Ref);

with qw(Cache::Ref::Role::API);

sub get { return }
sub hit { return }
sub set { return }
sub remove { return }
sub clear { return }
sub expire { return }
sub compute { return }

__PACKAGE__->meta->make_immutable;

__PACKAGE__;



=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Null - Caches nothing

=head1 SYNOPSIS

    # useful for comparing the effect of a cache compared to no
    # caching without code changes:

    my $c = Cache::Profile::Compare->new(
        caches => [
            Cache::Ref::Foo->new( ... ),
            Cache::Ref->Null->new,
        ],
    );

=head1 DESCRIPTION

This cache implementation will cache nothing.

This is primarily intended for testing or comparing runtime
without a cache against runtime with a cache.

It's like L<Cache::Null> but supports the additional methods in
L<Cache::Ref>.

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut


__END__


# ex: set sw=4 et:
