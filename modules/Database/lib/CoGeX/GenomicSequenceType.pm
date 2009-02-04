package CoGeX::GenomicSequenceType;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
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
__PACKAGE__->has_many("dataset_groups"=>"CoGeX::DatasetGroup","genomic_sequence_type_id");


sub resolve : ResultSet {
    my $self = shift;
    my $info = shift;
    return $info if ref($info) =~ /GenomicSequenceType/;
    return $self->find($info) if $info =~ /^\d+$/;
    my @res = $self->search({
			     'name' => { '-like' => $info . '%'}, 
			    }
			    ,{});
    @res = $self->search({
			     'name' => { '-like' => '%' . $info . '%'}, 
			    }
			    ,{}) unless scalar @res;
    return wantarray? @res : \@res;
}


1;

