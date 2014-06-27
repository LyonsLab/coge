package CoGeX::Result::AnnotationType;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::Result::AnnotationType

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<annotation_type> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<annotation_type_id> B<Primary Key>
Type:INT, Default: undef, Nullable: no, Size: 11

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255

C<annotation_type_group_id>
Type: INT, Default: undef, Nullable: yes, Size: 11

Belongs to CCoGeX::Result::AnnotationTypeGroup> via C<annotation_type_group_id>
Has many CCoGeX::Result::Annotation> via C<annotation_type_id>
Has many CCoGeX::Result::ExperimentAnnotation> via C<annotation_type_id>

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("annotation_type");
__PACKAGE__->add_columns(
  "annotation_type_id",
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
  "annotation_type_group_id",
  { data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
);
__PACKAGE__->set_primary_key("annotation_type_id");

__PACKAGE__->has_many("feature_annotations" => "CoGeX::Result::FeatureAnnotation", 'annotation_type_id');
__PACKAGE__->has_many("experiment_annotations" => "CoGeX::Result::ExperimentAnnotation", 'annotation_type_id');
__PACKAGE__->has_many("list_annotations" => "CoGeX::Result::ListAnnotation", 'annotation_type_id');

__PACKAGE__->belongs_to("annotation_type_group" => "CoGeX::Result::AnnotationTypeGroup", 'annotation_type_group_id');

################################################ subroutine header begin ##

=head2 group

 Usage     : $AnnotationType_obj->group->AnnotationTypeGroup_object_method_or_value
 Purpose   : Returns the AnnotationTypeGroup object associated with this AnnotationType object.
 Returns   : AnnotationTypeGroup object
 Argument  :
 Throws    : None
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub group
  {
    shift->annotation_type_group(@_);
  }

sub desc
    {
      shift->description(@_);
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
