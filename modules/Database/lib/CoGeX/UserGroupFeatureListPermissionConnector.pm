package CoGeX::UserGroupFeatureListPermissionConnector;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("user_group_feature_list_permission_connector");
__PACKAGE__->add_columns(
  "user_group_feature_list_permission_connector_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_group_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "feature_list_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "permission_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
);
__PACKAGE__->set_primary_key("user_group_feature_list_permission_connector_id");

1;

