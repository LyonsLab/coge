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
    return $self->find($info) if $info =~ /^\d+$/;
    return $self->search({
			  'name' => { '-like' => '%' . $info . '%'}, 
			 }
			 ,{});
}


sub current_datasets
  {
    my $self = shift;
    my %opts = @_;
    my $type = $opts{type} || $opts{genomic_sequence_type} || $opts{sequence_type};
    $type =1 unless $type;
    my %data;
    my $typeid;
    $typeid = ref($type) =~/Type/ ? $type->id : $type;
    my $version;
    ds_loop: foreach my $ds ($self->datasets({},{distict=>'version',order_by=>'version desc'}))
      {
	next unless $ds->sequence_type && $ds->sequence_type->id eq $typeid;
	$version = $ds->version unless $version;
	$version = $ds->version if $ds->version > $version;
#	next unless $version == $ds->version;
	my @chrs = $ds->get_chromosomes;
	foreach my $chr (@chrs)
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

sub genomic_sequence_types
  {
    my $self = shift;
    my %data;
    foreach my $ds ($self->datasets)
      {
	my $type = $ds->sequence_type;
	$data{$type->id} = $type;
      }
    return wantarray ? values %data : [values %data];
  }

1;
