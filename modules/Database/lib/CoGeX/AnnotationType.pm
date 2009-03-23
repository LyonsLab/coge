package CoGeX::AnnotationType;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::AnnotationType

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


Belongs to C<CoGeX::AnnotationTypeGroup> via C<annotation_type_group_id>
Has many C<CoGeX::Annotation> via C<annotation_type_id>

=head1 USAGE

 use CoGeX;
 
=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("annotation_type");
__PACKAGE__->add_columns(
  "annotation_type_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 100 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "annotation_type_group_id",
  { data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
);
__PACKAGE__->set_primary_key("annotation_type_id");

__PACKAGE__->has_many("annotations" => "CoGeX::Annotation", 'annotation_type_id');

__PACKAGE__->belongs_to("annotation_type_group" => "CoGeX::AnnotationTypeGroup", 'annotation_type_group_id');



################################################ subroutine header begin ##

=head2 group

 Usage     : $AnnotationType_obj->group->AnnotationTypeGroup_object_method_or_value
 Purpose   : Returns the AnnotationTypeGroup object associated with this AnnotationType object.
 Returns   : AnnotationTypeGroup object
 Argument  : None
 Throws    : None
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub group
  {
    shift->annotation_type_group(@_);
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
