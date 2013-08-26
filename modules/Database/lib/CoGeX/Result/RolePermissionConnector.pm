package CoGeX::Result::RolePermissionConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::RolePermissionConnector

=cut

__PACKAGE__->table("role_permission_connector");

=head1 ACCESSORS

=head2 role_permission_connector_id

  data_type: 'bigint'
  is_auto_increment: 1
  is_nullable: 0

=head2 role_id

  data_type: 'bigint'
  is_nullable: 1

=head2 permission_id

  data_type: 'bigint'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
    "role_permission_connector_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "role_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "permission_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
);
__PACKAGE__->set_primary_key("role_permission_connector_id");
__PACKAGE__->belongs_to( "role" => "CoGeX::Result::Role", "role_id" );
__PACKAGE__->belongs_to(
    "permission" => "CoGeX::Result::Permission",
    "permission_id"
);

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
