package CoGeX::GenomicSequenceType;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::

=head1 SYNOPSIS

  use CoGeX::
This object uses the DBIx::Class to define an interface to the C<genomic_sequence_type> table in the CoGe database.


=head1 DESCRIPTION


Has columns:
C<genomic_sequence_type_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255


Has many C<CoGeX::DatasetGroup> via C<genomic_sequence_type_id>

=head1 USAGE

=head1 METHODS

=cut


__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("genomic_sequence_type");
__PACKAGE__->add_columns(
  "genomic_sequence_type_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
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
__PACKAGE__->set_primary_key("genomic_sequence_type_id");
__PACKAGE__->has_many("dataset_groups"=>"CoGeX::DatasetGroup","genomic_sequence_type_id");


################################################ subroutine header begin ##

=head2 resolve

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

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


=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

=cut
