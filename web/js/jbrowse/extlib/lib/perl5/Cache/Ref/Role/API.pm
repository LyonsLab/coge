package Cache::Ref::Role::API;
BEGIN {
  $Cache::Ref::Role::API::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Role::API::VERSION = '0.04';
}
use Moose::Role;

use Carp qw(croak);

use namespace::autoclean;

requires qw(
    get
    set
    remove
    clear
    hit
    expire
);

sub compute {
    my ( $self, $key, $code ) = @_;

    croak "must specify key and code"
        unless defined($key) && defined($code);

    if ( defined( my $cached = $self->get($key) ) ) {
        return $cached;
    } else {
        my $value = $code->();
        $self->set( $key => $value );
        return $value;
    }
}

# ex: set sw=4 et:

__PACKAGE__;


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Role::API

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

