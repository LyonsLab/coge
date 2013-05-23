#!/usr/bin/perl

package Hash::Util::FieldHash::Compat;

use strict;
use warnings;

use constant REAL_FIELDHASH => do { local $@; eval { require Hash::Util::FieldHash } };

BEGIN {
	if ( REAL_FIELDHASH ) {
		require Hash::Util::FieldHash;
		Hash::Util::FieldHash->import(":all");
	} else {
		require Hash::Util::FieldHash::Compat::Heavy;
	}
}

our $VERSION = "0.03";

sub import {
	if ( REAL_FIELDHASH ) {
		my $class = "Hash::Util::FieldHash";

		shift @_;
		unshift @_, $class;

		goto $class->can("import");
	} else {
		my $class = shift;
		$class->export_to_level(1, $class, @_);
	}
}

__PACKAGE__

__END__

=pod

=head1 NAME

Hash::Util::FieldHash::Compat - Use L<Hash::Util::FieldHash> or ties, depending
on availability.

=head1 SYNOPSIS

	use Hash::Util::FieldHash::Compat;

	# pretend you are using L<Hash::Util::FieldHash>
	# under older perls it'll be Tie::RefHash::Weak instead (slower, but same behavior)

=head1 DESCRIPTION

Under older perls this module provides a drop in compatible api to
L<Hash::Util::FieldHash> using L<perltie>. When L<Hash::Util::FieldHash> is
available it will use that instead.

This way code requiring field hashes can benefit from fast, robust field hashes
on Perl 5.10 and newer, but still run on older perls that don't ship with that
module.

See L<Hash::Util::FieldHash> for all the details of the API.

=head1 SEE ALSO

L<Hash::Util::FieldHash>, L<Tie::RefHash>, L<Tie::RefHash::Weak>.

=head1 VERSION CONTROL

This module is maintained using Darcs. You can get the latest version from
L<http://nothingmuch.woobling.org/code>, and use C<darcs send> to commit
changes.

=head1 AUTHOR

Yuval Kogman E<lt>nothingmuch@woobling.orgE<gt>

=head1 COPYRIGHT

	Copyright (c) 2008 Yuval Kogman. All rights reserved
	This program is free software; you can redistribute
	it and/or modify it under the same terms as Perl itself.

=cut
