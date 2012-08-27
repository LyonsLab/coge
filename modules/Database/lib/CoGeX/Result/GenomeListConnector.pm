package CoGeX::Result::GenomeListConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::Result::GenomeListConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<genome_list_connector> table in the CoGe database.
The C<genome_list_connector> table is used to associate C<genome> (genome) records with C<list> records.

=head1 DESCRIPTION

Has columns:
C<genome_list_connector_id> (Primary Key)
Type: INT, Default: yes, Nullable: no, Size: 11
Primary identification key for table.

C<genome_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<genome> table.

C<list_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<list> table.

Belongs to CCoGeX::Result::Genome> via C<genome_id>
Belongs to CCoGeX::Result::List> via C<list_id>

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("genome_list_connector");
__PACKAGE__->add_columns(
  "genome_list_connector_id",
  { data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
  "genome_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
  "list_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("genome_list_connector_id");

__PACKAGE__->belongs_to("genome" => "CoGeX::Result::Genome", "genome_id");
__PACKAGE__->belongs_to("list" => "CoGeX::Result::List", "list_id");

1;

=head1 BUGS

=head1 SUPPORT

=head1 AUTHORS

 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
