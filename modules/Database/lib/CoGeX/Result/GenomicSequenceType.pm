package CoGeX_dev::Result::GenomicSequenceType;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';
use lib '/home/mbomhoff/CoGe/Accessory/lib'; #FIXME 8/2/12 remove
use lib '/home/mbomhoff/CoGeX/lib'; #FIXME 8/2/12 remove
use CoGeX_dev::ResultSet::GenomicSequenceType;

=head1 NAME

CoGeX_dev::GenomicSequenceType

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<genomic_sequence_type> table in the CoGe database.
The C<genomic_sequence_type> table defines the name and description of a genomic sequence type, and is used by the C<dataset_group> table.
This table should not be confused with the C<sequence_type> table - that table defines sequence types for the C<sequence> table which is used by the C<feature> table.

=head1 DESCRIPTION

Has columns:
C<genomic_sequence_type_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100
Name of genomic sequence type.

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Description of the genomic sequence type.

Relates to CCoGeX_dev::Result::Genome> via C<genomic_sequence_type_id> in a one-to-many relationship.

=head1 USAGE

  use CoGeX_dev;

=head1 METHODS

=cut


__PACKAGE__->table("genomic_sequence_type");
__PACKAGE__->add_columns(
  "genomic_sequence_type_id",
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
__PACKAGE__->set_primary_key("genomic_sequence_type_id");
__PACKAGE__->has_many("genomes"=>"CoGeX_dev::Result::Genome","genomic_sequence_type_id");

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
