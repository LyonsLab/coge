package CoGeX::Result::Location;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::Location

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<location> table in the CoGe database.
The C<location> table describes the location of a feature (a record in the C<feature> table) on a give genomic sequence. This includes the beginning and end of the feature, the strand it exists on, and the name of the chromosome it exists on.

=head1 DESCRIPTION

Has columns:
C<location_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<start>
Type: INT, Default: 0, Nullable: no, Size: 11
Start location of feature (in terms of basepairs) on the genomic data.

C<stop>
Type: INT, Default: 0, Nullable: no, Size: 11
Stop location of feature (in terms of basepairs) on the genomic data.

C<chromosome>
Type: VARCHAR, Default: "", Nullable: no, Size: 255
Name of chromosome this location resides in.

C<feature_id>
Type: INT, Default: 0, Nullable: no, Size: 11
Reference to record in C<feature> table.

C<strand>
Type: TINYINT, Default: "", Nullable: no, Size: 4
Which strand the feature exists on, the options being Top (+) and Bottom (-).
[Editors note: Why is this a TINYINT?]

Belongs to CCoGeX::Result::Feature> via C<feature_id>

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("location");
__PACKAGE__->add_columns(
  "location_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "start",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "stop",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "chromosome",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "strand",
  { data_type => "TINYINT", default_value => "", is_nullable => 0, size => 4 },
);
__PACKAGE__->set_primary_key("location_id");

__PACKAGE__->belongs_to("feature" => "CoGeX::Result::Feature", 'feature_id');
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
