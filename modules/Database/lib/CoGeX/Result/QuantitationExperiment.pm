package CoGeX::Result::QuantitationExperiment;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::QuantitationExperiment

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<quantitation_experiment> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<quantitation_experiment_id> B<Primary Key>
Type:INT, Default: undef, Nullable: no, Size: 11

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255
Type: INT, Default: undef, Nullable: yes, Size: 11



Has many CCoGeX::Result::Quantitation> via C<quantitation_experiment_id>

=head1 USAGE

 use CoGeX;
 
=head1 METHODS

=cut

__PACKAGE__->table("quantitation_experiment");
__PACKAGE__->add_columns(
  "quantitation_experiment_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 1024 },
  "description",
  {
    data_type => "TEXT",
    default_value => undef,
    is_nullable => 1,
  },
);
__PACKAGE__->set_primary_key("quantitation_experiment_id");

__PACKAGE__->has_many("quantitations" => "CoGeX::Result::Quantitation", 'quantitation_experiment_id');



sub desc
    {
      shift->description(@_);
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
