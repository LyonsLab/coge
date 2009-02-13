package CoGeX::GenomicSequence;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("genomic_sequence");
__PACKAGE__->add_columns(
  "genomic_sequence_id",
  { data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
  "length",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
  "chromosome",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "dataset_group_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("genomic_sequence_id");
__PACKAGE__->belongs_to("dataset_group" => "CoGeX::DatasetGroup", "dataset_group_id");

1;

