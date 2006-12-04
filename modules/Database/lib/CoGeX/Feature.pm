package CoGeX::Feature;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("feature");
__PACKAGE__->add_columns(
  "feature_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "feature_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "dataset_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
);

__PACKAGE__->set_primary_key("feature_id");

# feature has many feature_names
__PACKAGE__->has_many( 'feature_names' => "CoGeX::FeatureName", 'feature_id');

# feature has many annotations
__PACKAGE__->has_many( 'annotations' => "CoGeX::Annotation", 'feature_id');

# feature has many locations
__PACKAGE__->has_many( 'locations' => "CoGeX::Location", 'feature_id');

# feature has many sequences
__PACKAGE__->has_many( 'sequences' => "CoGeX::Sequence", 'feature_id');


# feature_type has many features
__PACKAGE__->belongs_to("feature_type" => "CoGeX::FeatureType", 'feature_type_id');

# dataset has many features
__PACKAGE__->belongs_to("dataset" => "CoGeX::Dataset", 'dataset_id');

1;
