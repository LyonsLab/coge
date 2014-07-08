package CoGeX::Result::FeatureAnnotation;

use strict;
use warnings;

use base 'DBIx::Class::Core';
use CoGeX::ResultSet::FeatureAnnotation;

=head1 NAME

CoGeX::Annotation

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<annotation> table in the CoGe database.

=head1 DESCRIPTION

Object for interacting with the Annotations table in the CoGe database.
Has columns:

C<annotation_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<annotation>
Type: TEXT, Deault: "", Nullable: no, Size 65535

C<feature_id>
Type: INT, Default: 0, Nullable: no, Size: 11

C<annotation_type_id>
Type: INT, Default: 0, Nullable: no, Size: 11

Belongs to CCoGeX::Result::AnnotationType> via C<annotation_type_id>

Belongs to CCoGeX::Result::Feature> via C<feature_id>

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("feature_annotation");
__PACKAGE__->add_columns(
  "feature_annotation_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "annotation",
  { data_type => "TEXT", default_value => "", is_nullable => 0, size => 65535 },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "annotation_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "link",
  { data_type => "VARCHAR", default_value => "", is_nullable => 1, size => 1024 },
);
__PACKAGE__->set_primary_key("feature_annotation_id");

__PACKAGE__->belongs_to( annotation_type => 'CoGeX::Result::AnnotationType', 'annotation_type_id');
__PACKAGE__->belongs_to( feature => 'CoGeX::Result::Feature', 'feature_id');
__PACKAGE__->has_one("feature_name" => "CoGeX::Result::FeatureName", {'foreign.feature_id'=>'self.feature_id'});

################################################ subroutine header begin ##

=head2 type

 Usage     : $Annotation_obj->type->AnnotationType_object_method_or_value
 Purpose   : Shorthand for getting an Annotation Type from an Annotation object.
 Returns   : AnnotationType object.
 Argument  :
 Throws    : None.
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub type
  {
    shift->annotation_type(@_);
  }

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
