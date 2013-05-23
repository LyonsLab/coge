package Heap::Simple::XS;
use strict;
# use warnings;
use vars qw($VERSION);

$VERSION = "0.10";

require XSLoader;
XSLoader::load('Heap::Simple::XS', $VERSION);

sub implementation() {
    return __PACKAGE__;
}

1;
__END__

=head1 NAME

Heap::Simple::XS - An XS implementation of the Heap::Simple interface

=head1 SYNOPSIS

    # Let Heap::Simple decide which implementation that provides its interface
    # it will load and use. This may be Heap::Simple::XS or it may not be.
    # Still, this is the normal way of using Heap::Simple
    use Heap::Simple;
    my $heap = Heap::Simple->new(...);
    # Use heap as described in the Heap::Simple documentation

    # If for some reason you insist on using this version:
    use Heap::Simple::XS;
    my $heap = Heap::Simple::XS->new(...);
    # Use the XS heap as described in the Heap::Simple documentation

=head1 DESCRIPTION

This module provides an XS implementation of the interface described
in L<Heap::Simple|Heap::Simple>. Look there for a description.

=head1 NOTES

=over

=item

Even though this implementation is written in C, it fully supports
overloading and magic (like L<ties|perltie>).

=item

The dirty option will do several things.

=over 4

=item

It will cause scalars for the C<E<lt>> and C<E<gt>> orders
to be stored internally as an NV (double or long double). This means you lose
magic, overload and any internal integer representation.

=item

The C<E<lt>> and C<E<gt>> order will cause C<Array> and C<Hash> elements
to get their key internally cached as an NV. So indirect changes to the value
won't be noticed anymore (but most of the time you shouldn't do that anyways).
It also means these will start behaving like a wrapped heap type, so they
return true for L<wrapped|Heap::Simple/wrapped> and support
L<key_insert|Heap::Simple/key_insert> and
L<key_absorb|Heap::Simple/key_absorb>.

=item

The C<E<lt>> and C<E<gt>> order will cause C<Object> and C<Any> elements
to store the key as an NV (these two already were wrapped heap types).

=item

It has no effect on C<Method> and C<Function> element types since it's
assumed you B<want> the key recalculations for these for some reason (if you
didn't, you would have asked for C<Object> or C<Any> elements).

=back

=item

Heap::Simple->implementation will return C<"Heap::Simple::XS"> if it selected
this module.

=back

=head1 EXPORT

None.

=head1 SEE ALSO

L<Heap::Simple>,
L<Heap::Simple::Perl>

=head1 AUTHOR

Ton Hospel, E<lt>Heap-Simple-XS@ton.iguana.beE<gt>

Parts are inspired by code by Joseph N. Hall
L<http://www.perlfaq.com/faqs/id/196>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2004 by Ton Hospel

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.6.1 or,
at your option, any later version of Perl 5 you may have available.

=cut
