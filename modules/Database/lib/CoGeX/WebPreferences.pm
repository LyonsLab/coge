package CoGeX::WebPreferences;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("web_preferences");
__PACKAGE__->add_columns(
  "id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "page",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "options",
  { data_type => "TEXT", default_value => "", is_nullable => 0 },
);
__PACKAGE__->set_primary_key("user_id");

1;

