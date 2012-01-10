package CoGeX::Result::List;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::FeatureList

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<feature_list> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<feature_list_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 50

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255

C<feature_list_group_id>
Type: INT, Default: undef, Nullable: yes, Size: 10

C<notes>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 1024


=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("list");
__PACKAGE__->add_columns(
  "list_id",
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
  "list_group_id",
  { data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
  "notes",
  {
    data_type => "TEXT",
    default_value => undef,
    is_nullable => 1,
  },
  "link",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 1024,
  },
  "user_id",
  {
   data_type => "INT",
   default_value => undef,
   is_nullable => 0,
   size => 11
  },
  "public",
  {
   data_type => "BOOLEAN",
   default_value => 0,
   is_nullable => 0,
   size => 1
  },
);
__PACKAGE__->set_primary_key("list_id");

# feature_type has many features
__PACKAGE__->has_many("list_connectors" => "CoGeX::Result::ListConnector", 'list_id');

# dataset has many features
__PACKAGE__->belongs_to("list_group" => "CoGeX::Result::ListGroup", 'list_group_id');


sub features
  {
    my $self = shift;
    my %opts = @_;
    my @feats;
    foreach my $conn ($self->list_connectors)
      {
	push @feats, $conn->feature if $conn->feature_id;
      }
    return wantarray ? @feats : \@feats;
  }

sub genomes
  {
    my $self = shift;
    my %opts = @_;
    my @dsgs;
    foreach my $conn ($self->list_connectors)
      {
	push @dsgs, $conn->genome if $conn->dataset_group_id;
      }
    return wantarray ? @dsgs : \@dsgs;
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
