package CoGeX::Result::Role;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

My::Schema::Result::Role

=cut

__PACKAGE__->table("role");

=head1 ACCESSORS

=head2 role_id

  data_type: 'bigint'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 description

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "role_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "description",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);
__PACKAGE__->set_primary_key("role_id");
# Role has many permissions
__PACKAGE__->has_many('role_permission_connectors'=>"CoGeX::Result::RolePermissionConnector",'role_id');
__PACKAGE__->has_many('user_groups'=>'CoGeX::Result::UserGroup', 'role_id');


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2011-08-29 09:28:12
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:xlyHeSUJRzYisE2OZEwHTQ


# You can replace this text with custom content, and it will be preserved on regeneration
1;
