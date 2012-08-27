package CoGeX::Result::ListCollection;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::Result::ListCollection

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<feature_list_group> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<feature_list_group_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 50

C<description>
Type: VARCHAR, Default: "", Nullable: no, Size: 255

C<notes>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 1024

C<link>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 1024

C<restricted>
Type: BOOLEAN, Default: 0, Nullable: no, Size: 1

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("list_collection");
__PACKAGE__->add_columns(
  "list_collection_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "description",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 1024 },
  "link",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 1024,
  },
  "restricted",
  {
   data_type => "BOOLEAN",
   default_value => 0,
   is_nullable => 0,
   size => 1
  },

);
__PACKAGE__->set_primary_key("list_collection_id");
__PACKAGE__->has_many("list_collection_connectors" => "CoGeX::Result::ListCollectionConnector", 'list_collection_id');

sub lists
  {
    my $self = shift;
    my %opts = @_;
    my @lists;
    foreach my $conn ($self->list_collection_connectors)
      {
	push @lists, $conn->list;
      }
    return wantarray ? @lists : \@lists;

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
