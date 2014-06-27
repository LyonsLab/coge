package CoGeX::Result::WebPreferences;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::WebPreferences

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<web_preferences> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<id>
Type: INT, Default: undef, Nullable: no, Size: 10

C<user_id> (Primary Key)
Type: INT, Default: "", Nullable: no, Size: 10

C<page>
Type: VARCHAR, Default: "", Nullable: no, Size: 255

C<options>
Type: TEXT, Default: "", Nullable: no, Size: N/A

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("web_preferences");
__PACKAGE__->add_columns(
  "id",  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "page",  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "options",  { data_type => "TEXT", default_value => "", is_nullable => 0 },
);
__PACKAGE__->set_primary_key("user_id");

1;

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
