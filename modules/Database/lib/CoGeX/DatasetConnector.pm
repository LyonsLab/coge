package CoGeX::DatasetConnector;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("dataset_connector");
__PACKAGE__->add_columns(
  "dataset_connector_id",
  { data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
  "dataset_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
  "dataset_group_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("dataset_connector_id");
__PACKAGE__->belongs_to("dataset_group" => "CoGeX::DatasetGroup", "dataset_group_id");
__PACKAGE__->belongs_to("dataset" => "CoGeX::Dataset", "dataset_id");

1;

