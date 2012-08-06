package CoGeX_dev::Result::Permission;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

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
__PACKAGE__->has_many('role_permission_connectors'=>"CoGeX_dev::Result::RolePermissionConnector","permission_id");

# Created by DBIx::Class::Schema::Loader v0.07002 @ 2011-08-29 09:28:12
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:+8kuH2wjE9fLJGNhIGz/Ww


# You can replace this text with custom content, and it will be preserved on regeneration
1;
