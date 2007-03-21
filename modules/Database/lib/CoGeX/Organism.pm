package CoGeX::Organism;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("organism");
__PACKAGE__->add_columns(
  "organism_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 100 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
);
__PACKAGE__->set_primary_key("organism_id");

__PACKAGE__->has_many("datasets" => "CoGeX::Dataset", 'organism_id');

sub resolve : ResultSet {
    my $self = shift;
    my $info = shift;
    return $info if ref($info) =~ /Organism/;
    return $self->search({
			  '-or'=>[
				  { 'name' => { '-like' => '%' . $info . '%'}}, 
				  { 'organism_id' => $info }
				 ],
			 }
			 ,{});
}


1;
