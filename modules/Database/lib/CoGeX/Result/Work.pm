package CoGeX::Result::Work;

use strict;
use warnings;
use base 'DBIx::Class';

=head1 NAME

CoGeX::Work

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<work> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<id>
Type: INT, Default: undef, Nullable: no, Size: 10

C<user_id> (Primary Key)
Type: INT, Default: "", Nullable: no, Size: 10

C<page>
Type: VARCHAR, Default: "", Nullable: no, Size: 100

C<parameter>
Type: TEXT, Default: "", Nullable: yes, Size: N/A

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("work");
__PACKAGE__->add_columns(
  "work_id",  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "page",  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 100 },
  "parameter", { data_type => "TEXT", default_value => "", is_nullable => 1 },
  "name", { data_type => "VARCHAR", default_value => "", is_nullable => 1, size=>256 },
  "description", { data_type => "VARCHAR", default_value => "", is_nullable => 1, size=>1024 },
  "note", { data_type => "TEXT", default_value => "", is_nullable => 1},
  "archive", { data_type => "SMALLINT", default_value => "0", is_nullable => 1, size=>1 },
  "image_id", { data_type => "INT", default_value => "", is_nullable => 1, size=>10 },
  "date", { data_type => "TIMESTAMP", default_value => "", is_nullable => 0 },
  "link", { data_type => "VARCHAR", default_value => "", is_nullable => 1, size=>1024 },
);
__PACKAGE__->set_primary_key("work_id");
__PACKAGE__->belongs_to('user'=>"CoGeX::Result::User","user_id");
__PACKAGE__->belongs_to('image'=>"CoGeX::Result::Image","image_id");
__PACKAGE__->has_many('work_orders'=>"CoGeX::Result::WorkOrder","work_id");
1;

=head1 AUTHORS

 Josh Kane
 Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
