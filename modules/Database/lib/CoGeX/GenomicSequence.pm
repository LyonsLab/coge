package CoGeX::GenomicSequence;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("genomic_sequence");
__PACKAGE__->add_columns(
  "genomic_sequence_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "start",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "stop",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "chromosome",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 15 },
  "sequence_data",
  {
    data_type => "LONGTEXT",
    default_value => "",
    is_nullable => 0,
    size => 4294967295,
  },
  "dataset_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "genomic_sequence_type_id",
  { data_type => "INT", default_value => 1, is_nullable => 0, size => 10 },
);
__PACKAGE__->set_primary_key("genomic_sequence_id");

__PACKAGE__->belongs_to("dataset" => "CoGeX::Dataset", "dataset_id");
__PACKAGE__->belongs_to("genomic_sequence_type" => "CoGeX::GenomicSequenceType","genomic_sequence_type_id");

1;

