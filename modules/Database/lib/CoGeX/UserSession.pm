package CoGeX::UserSession;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("user_session");
__PACKAGE__->add_columns(
  "user_session_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "date",
  { data_type => "DATETIME", default_value => "", is_nullable => 0, size => 19 },
);
__PACKAGE__->set_primary_key("user_session_id");




1;

