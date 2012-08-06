package CoGeX_dev::Result::ListCollectionImageConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX_dev::Result::ListCollectionImageConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<list_collection_image_connector> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<list_collection_image_connector_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<list_collection_id>
Type: INT, Default: "", Nullable: no, Size: 10

C<image_id>
Type: INT, Default: "", Nullable: no, Size: 10


=head1 USAGE

  use CoGeX_dev;

=head1 METHODS

=cut

__PACKAGE__->table("list_collection_image_connector");
__PACKAGE__->add_columns(
  "list_collection_image_connector_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "list_collection_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
  "image_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("list_collection_image_connector_id");
__PACKAGE__->belongs_to("list_collection" => "CoGeX_dev::Result::ListCollection", 'list_collection_id');
__PACKAGE__->belongs_to("image" => "CoGeX_dev::Result::Image", 'image_id');
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
