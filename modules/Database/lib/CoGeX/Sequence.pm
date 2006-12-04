package CoGeX::Sequence;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("sequence");
__PACKAGE__->add_columns(
  "sequence_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "sequence_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "sequence_data",
  {
    data_type => "LONGTEXT",
    default_value => "",
    is_nullable => 0,
    size => 4294967295,
  },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
);

__PACKAGE__->set_primary_key("sequence_id");

__PACKAGE__->belongs_to("feature" => "CoGeX::Feature", 'feature_id');
__PACKAGE__->belongs_to("sequence_type" => "CoGeX::SequenceType", 'sequence_type_id');
1;

