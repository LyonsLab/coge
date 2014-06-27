package CoGeX::Result::AnnotationTypeGroup;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::AnnotationTypeGroup

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<annotation_type_group> table in the CoGe database.
The C<annotation_type_group> table contains the name and description of all the Annotation Type Groups.

=head1 DESCRIPTION

Has columns:
C<annotation_type_group_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100
Name of Annotation Type Group.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Description of Annotation Type Group.

Relates to CCoGeX::Result::AnnotationType> via C<annotation_type_group_id>; has a one-to-many relationship.

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("annotation_type_group");
__PACKAGE__->add_columns(
  "annotation_type_group_id",
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
__PACKAGE__->set_primary_key("annotation_type_group_id");

__PACKAGE__->has_many("annotation_types" => "CoGeX::Result::AnnotationType", 'annotation_type_group_id');

sub desc
  {
    shift->description(@_);
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
