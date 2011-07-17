package CoGeX::Result::DataSource;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::DataSource

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<data_source> table in the CoGe database.
The C<data_source> table contains infomation on the source of the data described by a record in the C<dataset> table.
This includes the name and description of the source, as well as URL to the origonal source.

=head1 DESCRIPTION


Has columns:
C<data_source_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100
Name of source.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Description of source.

C<link>
Type: TEXT, Default: undef, Nullable: yes, Size: 65535
URL to origonal source.

Relates to CCoGeX::Result::Dataset> via C<data_source_id>, one-to-many relationship.

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("data_source");
__PACKAGE__->add_columns(
  "data_source_id",
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
  "link",
  {
    data_type => "TEXT",
    default_value => undef,
    is_nullable => 1,
    size => 65535,
  },
);
__PACKAGE__->set_primary_key("data_source_id");


__PACKAGE__->has_many('datasets'=>"CoGeX::Result::Dataset", "data_source_id");

sub desc
  {
    shift->description(@_);
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
