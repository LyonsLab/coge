package CoGeX::FeatureListGroupImageConnector;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("feature_list_group_image_connector");
__PACKAGE__->add_columns(
  "feature_list_group_image_connector_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "feature_list_group_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "image_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
);
__PACKAGE__->set_primary_key("feature_list_group_image_connector_id");

1;

