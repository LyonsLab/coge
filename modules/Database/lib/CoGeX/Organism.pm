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


sub current_datasets
  {
    my $self = shift;
    my %data;
    my $version;
    ds_loop: foreach my $ds ($self->datasets({},{distict=>'version',order_by=>'version desc'}))
      {
	$version = $ds->version unless $version;
	foreach my $chr ($ds->get_chromosomes)
	  {
	    #this is a hack but the general problem is that some organisms have different chromosomes at different versions, however, partially complete genomes will have many contigs and different versions will have different contigs.  So, to get around this, there is a check to see if the chromosome name has contig in it, if so, then only the most current version is used.  Otherwise, all versions are game.
	    next unless $chr;
	    if ($chr =~ /contig/i)
	      {
		next ds_loop if $ds->version ne $version;
	      }
	    $data{$chr} = $ds unless $data{$chr};
	    $data{$chr} = $ds if $ds->version > $data{$chr}->version;
	  }
      }
    %data = map {$_->id,$_} values %data;
    return wantarray ? values %data : [values %data];
  }

1;
