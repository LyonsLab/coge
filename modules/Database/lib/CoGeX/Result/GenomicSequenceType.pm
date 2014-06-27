package CoGeX::Result::GenomicSequenceType;

use strict;
use warnings;

use base 'DBIx::Class::Core';
use CoGeX::ResultSet::GenomicSequenceType;
use JSON qw(encode_json);

=head1 NAME

CoGeX::GenomicSequenceType

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<genomic_sequence_type> table in the CoGe database.
The C<genomic_sequence_type> table defines the name and description of a genomic sequence type, and is used by the C<dataset_group> table.
This table should not be confused with the C<sequence_type> table - that table defines sequence types for the C<sequence> table which is used by the C<feature> table.

=head1 DESCRIPTION

Has columns:
C<genomic_sequence_type_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100
Name of genomic sequence type.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Description of the genomic sequence type.

Relates to CCoGeX::Result::Genome> via C<genomic_sequence_type_id> in a one-to-many relationship.

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("genomic_sequence_type");
__PACKAGE__->add_columns(
  "genomic_sequence_type_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 256 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 1024,
  },
);
__PACKAGE__->set_primary_key("genomic_sequence_type_id");
__PACKAGE__->has_many("genomes"=>"CoGeX::Result::Genome","genomic_sequence_type_id");

################################################ subroutine header begin ##

=head2 to_*

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub to_hash {
    my $self = shift;
    return {
        id => $self->id,
        name => $self->name,
        description => $self->description
    }
}

sub to_json {
    return encode_json( shift->to_hash );
}

################################################ subroutine header begin ##

=head2 info

 Usage     :
 Purpose   : provides quick information about the type
 Returns   : a string
 Argument  :
 Throws    :
 Comments  : name, description
           :

See Also   :

=cut

################################################## subroutine header end ##

sub info
{
	my $self = shift;
	my $info = $self->name;
	$info .= ": " . $self->description      if $self->description;
	return $info;
}

1;

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
