package CoGeX::Result::Image;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::Image

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<image> table in the CoGe database.

=head1 DESCRIPTION


Has columns:
C<image_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 50
Name of image.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 2048
Description of image.

C<image>
Type: LONGBLOB, Default: "", Nullable: no, Size: 4294967295
Blob containing image data.

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut


__PACKAGE__->table("image");
__PACKAGE__->add_columns(
  "image_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "filename",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size =>256 },
  "image",
  {
    data_type => "LONGBLOB",
    default_value => "",
    is_nullable => 0,
    size => 4294967295,
  },
);
__PACKAGE__->set_primary_key("image_id");
__PACKAGE__->has_one( list_annotation => 'CoGeX::Result::ListAnnotation', 'image_id');

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
