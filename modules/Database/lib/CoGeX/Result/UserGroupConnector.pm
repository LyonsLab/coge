package CoGeX::Result::UserGroupConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::UserGroupConnector

=cut

__PACKAGE__->table("user_group_connector");

=head1 ACCESSORS

=head2 user_group_connector_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 user_id

  data_type: 'integer'
  is_nullable: 0

=head2 user_group_id

  data_type: 'integer'
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "user_group_connector_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "user_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "user_group_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("user_group_connector_id");
__PACKAGE__->belongs_to("user"=>"CoGeX::Result::User",'user_id');
__PACKAGE__->belongs_to("user_group"=>"CoGeX::Result::UserGroup",'user_group_id');

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
