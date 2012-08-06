package CoGeX_dev::Result::RolePermissionConnector;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

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
			 { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
			 "role_id",
			 { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
			 "permission_id",
			 { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
			);
__PACKAGE__->set_primary_key("role_permission_connector_id");
__PACKAGE__->belongs_to("role"=>"CoGeX_dev::Result::Role","role_id");
__PACKAGE__->belongs_to("permission"=>"CoGeX_dev::Result::Permission","permission_id");

# Created by DBIx::Class::Schema::Loader v0.07002 @ 2011-08-29 09:28:12
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:1Wzqk5xNmf/0MiICncYHEw


# You can replace this text with custom content, and it will be preserved on regeneration
1;
