package CoGeX::Result::ExperimentAnnotation;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';
use CoGeX::ResultSet::ExperimentAnnotation;

=head1 NAME

CoGeX::ExperimentAnnotation

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<experiment_annotation> table in the CoGe database.

=head1 DESCRIPTION

Object for interacting with the Experiment Annotations table in the CoGe database.
Has columns:

C<experiment_annotation_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<annotation>
Type: TEXT, Deault: "", Nullable: no

C<experiment_id>
Type: INT, Default: 0, Nullable: no, Size: 11

C<annotation_type_id>
Type: INT, Default: 0, Nullable: no, Size: 11


Belongs to CCoGeX::Result::AnnotationType> via C<annotation_type_id>

Belongs to CCoGeX::Result::Experiment> via C<experiment_id>

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("experiment_annotation");
__PACKAGE__->add_columns(
  "experiment_annotation_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "annotation",
  { data_type => "TEXT", default_value => "", is_nullable => 0 },
  "experiment_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "annotation_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("experiment_annotation_id");

__PACKAGE__->belongs_to( annotation_type => 'CoGeX::Result::AnnotationType', 'annotation_type_id');
__PACKAGE__->belongs_to( experiment => 'CoGeX::Result::Experiment', 'experiment_id');



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


=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons
 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO
