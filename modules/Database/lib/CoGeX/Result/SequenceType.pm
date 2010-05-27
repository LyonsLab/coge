package CoGeX::Result::SequenceType;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::SequenceType

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<sequence_type> table in the CoGe database.
The C<sequence_type> table contains names and descriptions of the various sequence types stored in the C<sequence> table.

=head1 DESCRIPTION

Has columns:
C<sequence_type_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100
Name of sequence type.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Description of sequence type.

Relates to CCoGeX::Result::Sequence> via C<sequence_type_id>; one-to-many relationship.


=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut


__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("sequence_type");
__PACKAGE__->add_columns(
  "sequence_type_id",
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

__PACKAGE__->set_primary_key("sequence_type_id");

__PACKAGE__->has_many( 'sequences' => "CoGeX::Result::Sequence", 'sequence_type_id');

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
