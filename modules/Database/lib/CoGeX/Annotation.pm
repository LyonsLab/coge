package CoGeX::Annotation;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("annotation");
__PACKAGE__->add_columns(
  "annotation_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
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

sub esearch : ResultSet
  {
    my $self = shift;
    my $join = $_[1]{'join'};
    map { push(@$join, $_ ) }
        ('annotation_type');


    my $prefetch = $_[1]{'prefetch'};
    map { push(@$prefetch, $_ ) }
        ('annotation_type',
            { 'annotation_type' => 'annotation_type_group' }
        );

    $_[1]{'join'} = $join;
    $_[1]{'prefetch'} = $prefetch;
    return $self->search( @_ );

  }

sub type
  {
    shift->annotation_type(@_);
  }

1;
