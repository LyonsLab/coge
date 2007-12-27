package CoGeX::GenomicSequenceType;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("genomic_sequence_type");
__PACKAGE__->add_columns(
  "genomic_sequence_type_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 50 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
);
__PACKAGE__->set_primary_key("genomic_sequence_type_id");


__PACKAGE__->has_many("genomic_sequences"=>"CoGeX::GenomicSequence","genomic_sequence_type_id");

1;

