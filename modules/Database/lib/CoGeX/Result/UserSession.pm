package CoGeX::Result::UserSession;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::UserSession

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user_session> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<user_session_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<user_id>
Type: INT, Default: "", Nullable: no, Size: 10

C<date>
Type: DATETIME, Default: "", Nullable: no, Size: 19

C<session>
Type: VARCHAR, Default: "", Nullable: no, Size: 22

Belongs to CCoGeX::Result::User> via C<user_id>

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("user_session");
__PACKAGE__->add_columns(
  "user_session_id",#  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",#  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "date",#  { data_type => "DATETIME", default_value => "", is_nullable => 0},
  "session",#  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 22 },
);
__PACKAGE__->set_primary_key("user_session_id");
__PACKAGE__->belongs_to("user"=>"CoGeX::Result::User", 'user_id');

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
