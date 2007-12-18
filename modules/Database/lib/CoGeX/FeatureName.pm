package CoGeX::FeatureName;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("feature_name");
__PACKAGE__->add_columns(
  "feature_name_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 50 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "primary_name",
  { data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 1 },
);
__PACKAGE__->set_primary_key("feature_name_id");
__PACKAGE__->belongs_to("feature" => "CoGeX::Feature", "feature_id");


sub esearch : ResultSet {
    my $self = shift;
    my $join = $_[1]{'join'};

    my $prefetch = $_[1]{'prefetch'};
    map { push(@$prefetch, $_ ) } 
        ({ 'feature' => ['locations','feature_type'] });


    $_[1]{'join'} = $join;
    $_[1]{'prefetch'} = $prefetch;
    my $rs = $self->search(
         @_
    );
    return $rs;

}

1;

