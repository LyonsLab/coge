package CoGeX_dev::Result::ExperimentListConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX_dev::Result::ExperimentListConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<experiment_list_connector> table in the CoGe database.
The C<experiment_list_connector> table is used to associate C<experiment> records with C<list> records.

=head1 DESCRIPTION

Has columns:
C<experiment_list_connector_id> (Primary Key)
Type: INT, Default: yes, Nullable: no, Size: 11
Primary identification key for table.

C<experiment_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<experiment> table.

C<list_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<list> table.

Belongs to CCoGeX_dev::Result::Experiment> via C<experiment_id>
Belongs to CCoGeX_dev::Result::List> via C<list_id>

=head1 USAGE

 use CoGeX_dev;

=head1 METHODS

=cut

__PACKAGE__->table("experiment_list_connector");
__PACKAGE__->add_columns(
  "experiment_list_connector_id",
  { data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
  "experiment_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
  "list_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("experiment_list_connector_id");

__PACKAGE__->belongs_to("experiment" => "CoGeX_dev::Result::Experiment", "experiment_id");
__PACKAGE__->belongs_to("list" => "CoGeX_dev::Result::List", "list_id");

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
