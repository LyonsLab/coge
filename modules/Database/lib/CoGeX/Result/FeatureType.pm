package CoGeX::Result::FeatureType;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::FeatureType

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<feature_type> table in the CoGe database.
The C<feature_type> table describes a feature type. Pretty straight forward really.

=head1 DESCRIPTION

Has columns:
C<feature_type_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100
Name of the feature type.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Description of the feature type.

Relates to CCoGeX::Result::Feature> via C<feature_type_id> in a one-to-many relationship.

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("feature_type");
__PACKAGE__->add_columns(
  "feature_type_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 1024,
  },
);
__PACKAGE__->set_primary_key("feature_type_id");
__PACKAGE__->has_many("features"=>"CoGeX::Result::Feature","feature_type_id");

sub desc
{
  return shift->description(@_);
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
