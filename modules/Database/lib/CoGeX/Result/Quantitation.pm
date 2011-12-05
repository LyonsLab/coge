package CoGeX::Result::Quantitation;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';
#use CoGeX::ResultSet::Quantitation;

=head1 NAME

CoGeX::Quantitation

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<quantitation> table in the CoGe database.

=head1 DESCRIPTION

Object for interacting with the Quantitations table in the CoGe database.
Has columns:

C<quantitation_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<quantitation>
Type: TEXT, Deault: "", Nullable: no, Size 65535

C<feature_id>
Type: INT, Default: 0, Nullable: no, Size: 11

C<quantitation_type_id>
Type: INT, Default: 0, Nullable: no, Size: 11


Belongs to CCoGeX::Result::QuantitationType> via C<quantitation_type_id>

Belongs to CCoGeX::Result::Feature> via C<feature_id>

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("quantitation");
__PACKAGE__->add_columns(
  "quantitation_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "quantitation",
  { data_type => "FLOAT", default_value => "", is_nullable => 0 },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "quantitation_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("quantitation_id");

__PACKAGE__->belongs_to( quantitation_type => 'CoGeX::Result::QuantitationType', 'quantitation_type_id');
__PACKAGE__->belongs_to( feature => 'CoGeX::Result::Feature', 'feature_id');




################################################ subroutine header begin ##

=head2 type

 Usage     : $Quantitation_obj->type->QuantitationType_object_method_or_value
 Purpose   : Shorthand for getting an Quantitation Type from an Quantitation object.
 Returns   : QuantitationType object.
 Argument  : 
 Throws    : None.
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub type
  {
    shift->quantitation_type(@_);
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
