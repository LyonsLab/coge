package CoGeX::Result::DatasetConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<dataset_connector> table in the CoGe database.
The C<dataset_connector> table is used to associate C<genome> records with C<dataset> records.

=head1 DESCRIPTION

Has columns:
C<dataset_connector_id> (Primary Key)
Type: INT, Default: yes, Nullable: no, Size: 11
Primary identification key for table.

C<dataset_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<dataset> table.

C<genome_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<genome> table.

Belongs to CCoGeX::Result::Genome> via C<genome_id>
Belongs to CCoGeX::Result::Dataset> via C<dataset_id>

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("dataset_connector");
__PACKAGE__->add_columns(
  "dataset_connector_id",
  { data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
  "dataset_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
  "genome_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("dataset_connector_id");
__PACKAGE__->belongs_to("genome" => "CoGeX::Result::Genome", "genome_id");
__PACKAGE__->belongs_to("dataset" => "CoGeX::Result::Dataset", "dataset_id");

1;

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
