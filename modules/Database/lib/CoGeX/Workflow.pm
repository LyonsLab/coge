package CoGeX::Workflow;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::Workflow

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<web_preferences> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<id>
Type: INT, Default: undef, Nullable: no, Size: 10

C<user_id> (Primary Key)
Type: INT, Default: "", Nullable: no, Size: 10

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100

C<description>
Type: VARCHAR, Default: "", Nullable: yes, Size: 255

C<date>
Type: TIMESTAMP, Default: "time created", Nullable: no,


=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("workflow");
__PACKAGE__->add_columns(
  "workflow_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",{ data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "name", { data_type => "VARCHAR", default_value => "", is_nullable => 0, size=>100 },
  "description", { data_type => "VARCHAR", default_value => "", is_nullable => 1, size=>255 },
  "date", { data_type => "TIMESTAMP", default_value => "", is_nullable => 0 },
);
__PACKAGE__->set_primary_key("workflow_id");
__PACKAGE__->has_many('work_orders'=>"CoGeX::WorkOrder","work_id");
__PACKAGE__->belongs_to('user'=>"CoGeX::User","user_id");



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
