package CoGeX::Annotation;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("annotation");
__PACKAGE__->add_columns(
  "annotation_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "annotation",
  { data_type => "TEXT", default_value => "", is_nullable => 0, size => 65535 },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "annotation_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("annotation_id");

__PACKAGE__->belongs_to( annotation_type => 'CoGeX::AnnotationType', 'annotation_type_id');
__PACKAGE__->belongs_to( feature => 'CoGeX::Feature', 'feature_id');
1;

