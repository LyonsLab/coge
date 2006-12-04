package CoGeX::FeatureListConnector;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("feature_list_connector");
__PACKAGE__->add_columns(
  "feature_list_connector_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "feature_list_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "feature_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "preferred_name",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 50,
  },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 1024,
  },
);
__PACKAGE__->set_primary_key("feature_list_connector_id");

1;

