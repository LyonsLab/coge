package CoGeX::Location;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("location");
__PACKAGE__->add_columns(
  "location_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "start",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "stop",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "strand",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 2 },
  "chromosome",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 15 },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
);
__PACKAGE__->set_primary_key("location_id");

__PACKAGE__->belongs_to("feature" => "CoGeX::Feature", 'feature_id');
1;

