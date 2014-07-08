package CoGeX::Result::WorkOrder;

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::WorkOrder

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<work_order> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<id>
Type: INT, Default: undef, Nullable: no, Size: 10

C<workflow_id>
Type: INT, Default: "", Nullable: no, Size: 10

C<work_id>
Type: INT, Default: "", Nullable: no, Size: 10

C<work_order>
Type: INT, Default: "", Nullable: no, Size: 4

C<date>
Type: TIMESTAMP, Default: "time created", Nullable: no,

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("work_order");
__PACKAGE__->add_columns(
  "work_order_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "workflow_id",{ data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "work_id",{ data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "work_order", { data_type => "INT", default_value => "1", is_nullable => 0, size=>4 },
  "date", { data_type => "TIMESTAMP", default_value => "", is_nullable => 0 },
);
__PACKAGE__->set_primary_key("work_order_id");
__PACKAGE__->belongs_to('work'=>"CoGeX::Result::Work","work_id");
__PACKAGE__->belongs_to('workflow'=>"CoGeX::Result::Workflow","workflow_id");

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
