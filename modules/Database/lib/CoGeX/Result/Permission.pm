package CoGeX::Result::Permission;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::Permission

=cut

__PACKAGE__->table("permission");

=head1 ACCESSORS

=head2 permission_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 256

=head2 description

  data_type: 'varchar'
  is_nullable: 1
  size: 1024

=cut

__PACKAGE__->add_columns(
			 "permission_id",
			 { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
			 "name",
			 { data_type => "varchar", is_nullable => 0, size => 256 },
			 "description",
			 { data_type => "TEXT", is_nullable => 1, size => 1024 },
			);
__PACKAGE__->set_primary_key("permission_id");
__PACKAGE__->has_many('role_permission_connectors'=>"CoGeX::Result::RolePermissionConnector","permission_id");

1;

=head1 AUTHORS

 Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
