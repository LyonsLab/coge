package CoGe::Genome;
use strict;

use CoGe::Genome::DB::Annotation;
use CoGe::Genome::DB::Annotation_type;
use CoGe::Genome::DB::Annotation_type_group;
use CoGe::Genome::DB::Data_source;
use CoGe::Genome::DB::Data_information;
use CoGe::Genome::DB::Feature;
use CoGe::Genome::DB::Feature_name;
use CoGe::Genome::DB::Feature_type;
use CoGe::Genome::DB::Genomic_sequence;
use CoGe::Genome::DB::Location;
use CoGe::Genome::DB::Organism;
use CoGe::Genome::DB::Sequence;
use CoGe::Genome::DB::Sequence_type;
use CoGe::Genome::DB::User;
use CoGe::Genome::DB::User_session;
use CoGe::Genome::DB::User_group;
use CoGe::Genome::DB::User_group_connector;
use CoGe::Genome::DB::Permission;
use CoGe::Genome::DB::User_group_feature_list_permission_connector;
use CoGe::Genome::DB::Feature_list;
use CoGe::Genome::DB::Feature_list_group;
use CoGe::Genome::DB::Feature_list_connector;
use CoGe::Genome::DB::Feature_list_group_image_connector;
use CoGe::Genome::DB::Image;

use Carp qw(cluck);
use Data::Dumper;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();

}


########################################### main pod documentation begin ##



=head1 NAME

Genome - Genome

=head1 SYNOPSIS

  use Genome
  my $genome_obj = Genome->new();


=head1 DESCRIPTION

  Master object for accessing information in the Genome Database from the SBCS
  Group at UC Berkeley


=head1 USAGE
  


=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

############################################# main pod documentation end ##


################################################ subroutine header begin ##

=head2 new

 Usage     : my $genome_obj = CoGe::Genome->new();
 Purpose   : Creates a CoGe::Genome object
 Returns   : an CoGe::Genome object
 Argument  : none
 Throws    : none
 Comments  : creates a Genome object which CoGe will use to manipulate genome data,
	     annotations, and features

See Also   : CoGe::Genome

=cut

################################################## subroutine header end ##



sub new
{
    my ($class, %parameters) = @_;

    my $self = bless ({}, ref ($class) || $class);

    return ($self);
}

################################################ subroutine header begin ##

=head2 get_annotation_obj

 Usage     : $genome_obj->get_annotation_obj();
 Purpose   : gets an annotation object
 Returns   : an annotation object
 Argument  : none
 Throws    : 
 Comments  : creates an annotation object which has an annotation_id, an annotation,
 	     a feature_id, and an annotation_type_id

See Also   : CoGe::Genome::DB::Annotation

=cut

################################################## subroutine header end ##



sub get_annotation_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Annotation->new();
  }

################################################ subroutine header begin ##

=head2 get_anno_obj

 Usage     : $genome_obj->get_anno_obj();
 Purpose   : alias for get_annotation_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias

See Also   : get_annotation_obj

=cut

################################################## subroutine header end ##

sub get_anno_obj
  {
    my $self = shift;
    return $self->get_annotation_obj();
  }

################################################ subroutine header begin ##

=head2 get_annotation_type_obj

 Usage     : $genome_obj->get_annotation_type_obj();
 Purpose   : gets an annotation_type object
 Returns   : an annotation_type object
 Argument  : none
 Throws    : 
 Comments  : Creates an annotatio_type object which has an annotation_type_id,
 	     a name, a description, and an annotation_type_group_id

See Also   : CoGe::Genome::DB::Annotation_type

=cut

################################################## subroutine header end ##

sub get_annotation_type_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Annotation_type->new();
  }

################################################ subroutine header begin ##

=head2 get_anno_type_obj

 Usage     : $genome_obj->get_anno_type_obj();
 Purpose   : alias for get_annotation_type_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias

See Also   : get_annotation_type_obj

=cut

################################################## subroutine header end ##

sub get_anno_type_obj
  {
    my $self = shift;
    return $self->get_annotation_type_obj();
  }

################################################ subroutine header begin ##

=head2 get_annotation_type_group_obj

 Usage     : $genome_obj->get_annotation_type_group_obj();
 Purpose   : gets new Annotation_type_group object
 Returns   : an Annotation_type_group object
 Argument  : none
 Throws    : 
 Comments  : Creates an Annotation_type_group object which has an annotation_type_group_id,
 	     a name, and a description

See Also   : CoGe::Genome::DB::Annotation_type_group

=cut

################################################## subroutine header end ##

sub get_annotation_type_group_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Annotation_type_group->new();
  }

################################################ subroutine header begin ##

=head2 get_anno_type_obj

 Usage     : $genome_obj->get_anno_type_obj();
 Purpose   : alias for get_annotation_type_group_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get annotation_type_group_obj 

=cut

################################################## subroutine header end ##

sub get_anno_type_group_obj
  {
    my $self = shift;
    return $self->get_annotation_type_group_obj();
  }

################################################ subroutine header begin ##

=head2 get_data_source_obj

 Usage     : $genome_obj->get_data_source_obj();
 Purpose   : gets a new Data_source object 
 Returns   : a Data_source object 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Data_source

=cut

################################################## subroutine header end ##

sub get_data_source_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Data_source->new();
  }

################################################ subroutine header begin ##

=head2 get_source_obj

 Usage     : $genome_obj->get_source_obj();
 Purpose   : alias for get_data_source_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_data_source_obj 

=cut

################################################## subroutine header end ##

sub get_source_obj
  {
    my $self = shift;
    return $self->get_data_source_obj();
  }

################################################ subroutine header begin ##

=head2 get_dataset_obj

 Usage     : $genome_obj->get_dataset_obj();
 Purpose   : gets a Dataset object 
 Returns   : a Dataset object
 Argument  : none 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Dataset

=cut

################################################## subroutine header end ##

sub get_dataset_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Dataset->new();
  }


################################################ subroutine header begin ##

=head2 get_information_obj

 Usage     : $genome_obj->get_information_obj();
 Purpose   : alias for get_dataset_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : Data_information package changed to Dataset. This subroutine
 	     was modified to point to Dataset package, and warn that 
	     "data_information" is obselete terminology, and any call to it
	     should be modified to call get_dataset_obj.

See Also   : get_data_information_obj

=cut

################################################## subroutine header end ##

sub get_data_information_obj
  {
    my $self = shift;
    print STDERR "This method is obselete. Please use get_dataset_obj";
    return $self->get_dataset_obj();
  }

################################################ subroutine header begin ##

=head2 get_data_info_obj

 Usage     : $genome_obj->get_data_info_obj();
 Purpose   : alias for get_data_information_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_data_information_obj

=cut

################################################## subroutine header end ##

sub get_data_info_obj
  {
    my $self = shift;
    return $self->get_data_information_obj();
  }

################################################ subroutine header begin ##

=head2 get_dataset_obj

 Usage     : $genome_obj->get_dataset_obj();
 Purpose   : alias for get_data_information_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_data_information_obj

=cut

################################################## subroutine header end ##

sub get_dataset_obj
  {
    my $self = shift;
    return $self->get_data_information_obj();
  }

################################################ subroutine header begin ##

=head2 get_feature_obj

 Usage     : $genome_obj->get_feature_obj();
 Purpose   : Creates new feature object
 Returns   : A feature object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##

sub get_feature_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature->new();
  }

################################################ subroutine header begin ##

=head2 get_feat_obj

 Usage     : $genome_obj->get_feat_obj();
 Purpose   : alias for get_feature_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_feature_obj

=cut

################################################## subroutine header end ##

sub get_feat_obj
  {
    my $self = shift;
    return $self->get_feature_obj();
  }

################################################ subroutine header begin ##

=head2 get_feature_name_obj

 Usage     : $genome_obj->get_feature_name_obj();
 Purpose   : creates a Feature_name object 
 Returns   : a Feature_name object
 Argument  : none 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_name

=cut

################################################## subroutine header end ##

sub get_feature_name_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_name->new();
  }

################################################ subroutine header begin ##

=head2 get_feat_name

 Usage     : $genome_obj->get_feat_name_obj();
 Purpose   : alias for get_feature_name_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_feature_name_obj 

=cut

################################################## subroutine header end ##

sub get_feat_name_obj
  {
    my $self = shift;
    return $self->get_feature_name_obj();
  }

################################################ subroutine header begin ##

=head2 get_feature_type_obj

 Usage     : $genome_obj->get_feature_type_obj();
 Purpose   : Creates a get_feature_type object 
 Returns   : a get_feature_type object 
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_type

=cut

################################################## subroutine header end ##

sub get_feature_type_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_type->new();
  }

################################################ subroutine header begin ##

=head2 get_feat_type_obj

 Usage     : $genome_obj->get_feat_type_obj();
 Purpose   : alias for get_feature_type_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_feature_type_obj

=cut

################################################## subroutine header end ##

sub get_feat_type_obj
  {
    my $self = shift;
    return $self->get_feature_type_obj();
  }

################################################ subroutine header begin ##

=head2 get_genomic_sequence_obj

 Usage     : $genome_obj->get_genomic_sequence_obj();
 Purpose   : Creates a new Genomic_sequence object 
 Returns   : a Genomic_sequence object
 Argument  : none 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Genomic_sequence

=cut

################################################## subroutine header end ##

sub get_genomic_sequence_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Genomic_sequence->new();
  }

################################################ subroutine header begin ##

=head2 get_genomic_seq_obj

 Usage     : $genome_obj->get_genomic_seq_obj();
 Purpose   : alias for get_genomic_sequence_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_genomic_sequence_obj

=cut

################################################## subroutine header end ##

sub get_genomic_seq_obj
  {
    my $self = shift;
    return $self->get_genomic_sequence_obj();
  }

################################################ subroutine header begin ##

=head2 get_location_obj

 Usage     : $genome_obj->get_location_obj();
 Purpose   : Creates a Location object 
 Returns   : a Location object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Location

=cut

################################################## subroutine header end ##

sub get_location_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Location->new();
  }

################################################ subroutine header begin ##

=head2 get_loc_obj

 Usage     : $genome_obj->get_loc_obj();
 Purpose   : alias for get_location_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_location_obj

=cut

################################################## subroutine header end ##

sub get_loc_obj
  {
    my $self = shift;
    return $self->get_location_obj();
  }

################################################ subroutine header begin ##

=head2 get_organism_obj

 Usage     : $genome_obj->get_organism_obj();
 Purpose   : gets an Organism object
 Returns   : an Organism object
 Argument  : none
 Throws    : none
 Comments  : 

See Also   : CoGe::Genome::DB::Organism

=cut

################################################## subroutine header end ##

sub get_organism_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Organism->new();
  }

################################################ subroutine header begin ##

=head2 get_org_obj

 Usage     : $genome_obj->get_org_obj();
 Purpose   : alias for get_organism_obj
 Returns   :  
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_organism_obj 

=cut

################################################## subroutine header end ##

sub get_org_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Organism->new();
  }

################################################ subroutine header begin ##

=head2 get_species_obj

 Usage     : $genome_obj->get_species_obj();
 Purpose   : alias for get_organism_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_organism_obj

=cut
  
################################################## subroutine header end ##

sub get_species_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Organism->new();
  }

################################################ subroutine header begin ##

=head2 

 Usage     : $genome_obj->sequence();
 Purpose   : Creates a Sequence object  
 Returns   : a Sequence object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Sequence

=cut

################################################## subroutine header end ##

sub get_sequence_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Sequence->new();
  }


################################################ subroutine header begin ##

=head2 get_seq_obj

 Usage     : $genome_obj->get_seq_obj();
 Purpose   : alias for get_sequence_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_sequence_obj 

=cut

################################################## subroutine header end ##

sub get_seq_obj
  {
    my $self = shift;
    return $self->get_sequence_obj();
  }


################################################ subroutine header begin ##

=head2 get_sequence_type_obj

 Usage     : $genome_obj->get_sequence_type_obj();
 Purpose   : Creates a Sequence_type object 
 Returns   : a Sequence_type object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Sequence_type

=cut

################################################## subroutine header end ##

sub get_sequence_type_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Sequence_type->new();
  }


################################################ subroutine header begin ##

=head2 get_seq_type_obj

 Usage     : $genome_obj->get_seq_type_obj();
 Purpose   : alias for get_sequence_type_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : This is an alias 

See Also   : get_sequence_type_obj 

=cut

################################################## subroutine header end ##

sub get_seq_type_obj
  {
    my $self = shift;
    return $self->get_sequence_type_obj();
  }



################################################ subroutine header begin ##

=head2 get_user_obj

 Usage     : $genome_obj->get_user_obj();
 Purpose   : Creates a User object 
 Returns   : a User object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::User

=cut

################################################## subroutine header end ##

sub get_user_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::User->new();
  }

################################################ subroutine header begin ##

=head2 get_user_group_obj

 Usage     : $genome_obj->get_user_group_obj();
 Purpose   : Creates a User_group object 
 Returns   : a User_group object
 Argument  : none 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::User_group

=cut

################################################## subroutine header end ##

sub get_user_group_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::User_group->new();
  }

################################################ subroutine header begin ##

=head2 get_user_group_connector_obj

 Usage     : $genome_obj->get_user_group_connector_obj();
 Purpose   : Creates a User_group_connector object 
 Returns   : a User_group_connector object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::User_group_connector

=cut

################################################## subroutine header end ##

sub get_user_group_connector_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::User_group_connector->new();
  }

################################################ subroutine header begin ##

=head2 get_permission_obj

 Usage     : $genome_obj->get_permission_obj();
 Purpose   : Creates a Permission object 
 Returns   : a Permission object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Permission

=cut

################################################## subroutine header end ##

sub get_permission_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Permission->new();
  }

################################################ subroutine header begin ##

=head2 get_user_group_feature_list_permission_connector_obj

 Usage     : $genome_obj->get_user_group_feature_list_permission_connector_obj();
 Purpose   : Creates a User_group_feature_list_permission_connector object 
 Returns   : a User_group_feature_list_permission_connector object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::User_group_feature_list_permission_connector

=cut

################################################## subroutine header end ##

sub get_user_group_feature_list_permission_connector_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::User_group_feature_list_permission_connector->new();
  }

################################################ subroutine header begin ##

=head2 get_feature_list_obj

 Usage     : $genome_obj->get_feature_list_obj();
 Purpose   : Creates a Feature_list object 
 Returns   : a Feature_list object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_list

=cut

################################################## subroutine header end ##

sub get_feature_list_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_list->new();
  }

################################################ subroutine header begin ##

=head2 get_feature_list_group_obj

 Usage     : $genome_obj->get_feature_list_group_obj();
 Purpose   : Creates a Feature_list_group object 
 Returns   : a Feature_list_group object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_list_group

=cut

################################################## subroutine header end ##

sub get_feature_list_group_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_list_group->new();
  }

################################################ subroutine header begin ##

=head2 get_feature_list_connector_obj

 Usage     : $genome_obj->get_feature_list_connector_obj();
 Purpose   : Creates a Feature_list_connector object 
 Returns   : a Feature_list_connector object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_list_connector

=cut

################################################## subroutine header end ##

sub get_feature_list_connector_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_list_connector->new();
  }


################################################ subroutine header begin ##

=head2 get_user_session_obj

 Usage     : $genome_obj->get_user_session_obj();
 Purpose   : Creates a User_session object 
 Returns   : a User_session object
 Argument  : none
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::User_session

=cut

################################################## subroutine header end ##



sub get_user_session_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::User_session->new();
  }

################################################ subroutine header begin ##

=head2 get_image_obj

 Usage     : $genome_obj->get_image_obj();
 Purpose   : Creates an Image object 
 Returns   : an Image object
 Argument  : none
 Throws	   :
 Comments  :

See Also   : CoGe::Genome::DB::Image

=cut

################################################## subroutine header end ##

sub get_image_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Image->new();
  }


################################################ subroutine header begin ##

=head2 get_feature_list_group_image_connector_obj

 Usage     : $genome_obj->get_feature_list_group_image_connector_obj();
 Purpose   : Creates a Feature_list_group_image_connector object
 Returns   : a Feature_list_group_image_connector object
 Argument  : none
 Throws	   :
 Comments  :

See Also   : CoGe::Genome::DB::Feature_list_group_image_connector

=cut

################################################## subroutine header end ##

sub get_feature_list_group_image_connector_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_list_group_image_connector->new();
  }





=head1 FUNCTIONAL METHODS


################################################ subroutine header begin ##

=head2 get_features_by_name_and_info_id

 Usage     : my @feats = $genome_obj->get_features_by_name_and_info_id(name=>$name, id=>$id)
 Purpose   : get feature objects based on a name and a data_information id
 Returns   : wantarray of CoGe::Genome::Feature objects
 Argument  : a hash of key/value pairs
               name=> a string that is the name of a feature to be fetched by $self->get_feature_by_name
               id=>   an integer that is the database id from the data_information table
                      the id can be obtained from a data_information object
 Throws    : none
 Comments  : 

See Also   : $self->get_features_by_name
             CoGe::Data_information for information on that object

=cut

################################################## subroutine header end ##

sub get_features_by_name_and_info_id
  {
    my $self = shift;
    my %opts = @_;
    my ($name) = $opts{name} || $opts{NAME} ;
    my ($id) = $opts{info_id} || $opts{id} || $opts{ID};
    my @f;
    foreach my $feat ($self->get_features_by_name($name) )
      {	
	push @f, $feat if $feat->data_information->id eq $id;
      }
    return wantarray ? @f : \@f;
  }

################################################ subroutine header begin ##

=head2 get_feature_by_name_and_info_id

 Usage     : $genome_obj->get_feature_by_name_and_info_id($name);
 Purpose   : alias for get_features_by_name_and_info_id

=cut

################################################## subroutine header end ##


sub get_feature_by_name_and_info_id
  {
    my $self = shift;
    return $self->get_features_by_name_and_info_id(@_);
  }


################################################ subroutine header begin ##

=head2 get_features_by_name_and_dataset_id

 Usage     : $genome_obj->get_features_by_name_and_dataset_id($name);
 Purpose   : alias for get_features_by_name_and_info_id

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_dataset_id
  {
    my $self = shift;
    return $self->get_features_by_name_and_info_id(@_);
  }

################################################ subroutine header begin ##

=head2 get_feature_by_name_and_dataset_id

 Usage     : $genome_obj->get_feature_by_name_and_dataset_id($name);
 Purpose   : alias for get_features_by_name_and_info_id

=cut

################################################## subroutine header end ##


sub get_feature_by_name_and_dataset_id
  {
    my $self = shift;
    return $self->get_features_by_name_and_info_id(@_);
  }


################################################ subroutine header begin ##

=head2 get_features_by_name_and_version

 Usage     : my @feats = $genome_obj->get_features_by_name_and_version(name=>$name, version=>$version)
 Purpose   : get feature objects based on a name and a data_information version number
 Returns   : wantarray of CoGe::Genome::Feature objects
 Argument  : a hash of key/value pairs
               name=> a string that is the name of a feature to be fetched by $self->get_feature_by_name
               version=>   an integer.  This is obtained from a data_information object.
 Throws    : none
 Comments  : 

See Also   : sub get_feature_by_name
             CoGe::Data_information for information on that object

=cut

################################################## subroutine header end ##

sub get_features_by_name_and_version
  {
    my $self = shift;
    my @feats = $self->get_feature_obj->get_features_by_name_and_version(@_);
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 get_feature_by_name

 Usage     : $genome_obj->get_feature_by_name($name);
 Purpose   : get feature objects using the name of the feature
 Returns   : an array of feature objects (or nothing)
 Argument  : string / the name of the feature(s)
 Throws    : none
 Comments  : same as:
             map {$_->feature} $self->get_feat_name_obj->search(name=>$name)

See Also   : 

=cut

################################################## subroutine header end ##


sub get_feature_by_name
  {
    my $self = shift;
    my $name = shift;
    my @feats = map {$_->feature} $self->get_feat_name_obj->search(name=>$name);
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 get_features_by_name

 Usage     : $genome_obj->get_features_by_name($name);
 Purpose   : alias for get_feature_by_name
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_features_by_name
  {
    my $self = shift;
    my $name = shift;
    return $self->get_feature_by_name($name);
  }

################################################ subroutine header begin ##

=head2 get_feat_by_name

 Usage     : $genome_obj->get_feat_by_name($name);
 Purpose   : alias for get_feature_by_name
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feat_by_name
  {
    my $self = shift;
    my $name = shift;
    return $self->get_feature_by_name($name);
  }

################################################ subroutine header begin ##

=head2 get_feats_by_name

 Usage     : $genome_obj->get_feats_by_name($name);
 Purpose   : alias for get_feature_by_name
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feats_by_name
  {
    my $self = shift;
    my $name = shift;
    return $self->get_feature_by_name($name);
  }


################################################ subroutine header begin ##

=head2 get_protein_seq_by_feat_name

 Usage     : my @protein_seqs = $genome_obj->get_protein_seq_by_feat_name(name=>$name, version=>$version);
 Purpose   : Some features have an associated protein sequence which is stored
             in the sequence table with a sequence type "protein".  This routine 
             finds all features by a specified name and limited by a data_information 
             version and returns the first sequence for that feature it finds of type
             protein.
 Returns   : array of strings (protein sequence)
 Argument  : name => the name you wish to searhc (e.g. At3g010110)
             version => (optional) skips any feature whose version (from data_information)
                        does not equal this specified version.  (e.g. 6)
             hash   => Returns a hash (or hash ref) with keys as name of the feature that
                       matches the submitted name so that version numbers can be used
                       e.g. At3g010110.1=>MAR.. .  
             search_like => flag for "search like".  adds "%" after the accn name for the search
 Throws    : undef
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub get_protein_seq_by_feat_name
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name};
    my $version = $opts{version} || $opts{ver};
    my $hash = $opts{hash};
    my $search_like = $opts{search_like};
    return unless $name;
    my %seqs;
    my @seqs;
    my @no = $search_like ? $self->get_feature_name_obj->search_like(name=>$name."%") : $self->get_feature_name_obj->search(name=>$name);
    
    foreach my $no (@no)
      {
	if ($version)
	  {
	    next unless $no->feat->data_info->version eq $version;
	  }
	
	foreach my $seq ($no->feat->sequences)
	  {
	    push @seqs, $seq->sequence_data if $seq->seq_type->name =~ /prot/i;
	    $seqs{$no->name}=$seq->sequence_data if $seq->seq_type->name =~ /prot/i;
	  }
      }
    if ($hash)
      {
	return wantarray ? %seqs : \%seqs;
      }
    if (@seqs == 1)
      {
	return $seqs[0];
      }
    my @sorted = sort {length($b) <=> length ($a)} @seqs;
    return wantarray ? @sorted : \@sorted;
  }

################################################ subroutine header begin ##

=head2 get_prot_seq_by_feat_name

 Usage     : my @protein_seqs = $genome_obj->get_prot_seq_by_feat_name(name=>$name, version=>$version);
 Purpose   : Alias for get_protein_seq_by_feat_name

See Also   : get_protein_seq_by_feat_name

=cut

################################################## subroutine header end ##

sub get_prot_seq_by_feat_name
  {
    my $self = shift;
    return $self->get_protein_seq_by_feat_name(@_);
  }

################################################ subroutine header begin ##

=head2 get_protein_sequence_for_feature

 Usage     : my @protein_seqs = $genome_obj->get_protein_seq_by_feat($feat);
 Purpose   : Some features have an associated protein sequence which is stored
             in the sequence table with a sequence type "protein" (e.g. CDS sequences).  
             This routine finds and returns those sequences.
 Returns   : array of strings (protein sequence) or array ref (based on wantarray)
 Argument  : GoGe::Genomes::DB::Feature object
 Throws    : undef
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub get_protein_sequence_for_feature
  {
    my $self = shift;
    my $feat = shift;
    return unless $feat;
    my @seqs;
    foreach my $seq ($feat->sequences)
      {
	push @seqs, $seq->sequence_data if $seq->seq_type->name =~ /prot/i;
      }
    my @sorted = sort {length($b) <=> length ($a)} @seqs;
    return wantarray ? @sorted : \@sorted;
  }

################################################ subroutine header begin ##

=head2 get_genomic_seq_by_feat_name_and_type_name

 Usage     : my @seqs = $genome_obj->get_genomic_seq_by_feat_name_and_type_name(name=>$name, 
                                                                                type=>$type, 
                                                                                version=>$version);
 Purpose   : Gets the genomic sequence of a feature by the name of the feature and 
             by the type of the feature
 Returns   : array of strings (genomic DNA sequence)
 Argument  : name => the name you wish to searhc (e.g. At3g010110)
             type => the name of the feature type (e.g. CDS)
             version => (optional) skips any feature whose version (from data_information)
                        does not equal this specified version.  (e.g. 6)
             infoid  => (optional) skips any feature whose data information id 
                        does not match
 Throws    : undef, prints error to STDERR
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_genomic_seq_by_feat_name_and_type_name
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name};
    my $version = $opts{version} || $opts{ver};
    my $type = $opts{type};
    my $info_id = $opts{infoid};
    my %locs;
    foreach my $no ($self->get_feature_name_obj->search(name=>$name) )
      {
	foreach my $feat ($no->feat)
	  {
	    if ($version)
	      {
		next unless $feat->data_info->version eq $version;
	      }
	    if ($info_id)
	      {
		next unless $feat->data_info->id eq $info_id;
	      }
	    next unless $feat->feat_type->name =~ /$type/;
	    push @{$locs{$feat->id}}, $feat->locations;
	  }
      }
#    @locs = $self->get_location_obj->combine_overlaps(@locs);
    my $gso = $self->get_genomic_seq_obj;
    my $rev_seq_flag = 0;
    my @seqs;
    foreach my $id (keys %locs)
      {
	my $seq;
	
	foreach my $loc (sort {$a->start <=> $b->start} @{$locs{$id}})
	  {
	    $seq .= $gso->get_sequence(start=>$loc->start,
				       stop =>$loc->stop,
				       chr  =>$loc->chr,
				       strand=>$loc->strand,
				       org_id=>$loc->feat->org->id(),
				       data_info_id=>$loc->feat->data_info->id(),
				      );
	    $rev_seq_flag = 1 if $loc->strand =~ "-";
	  }
	$seq = reverse $seq if $rev_seq_flag;
	push @seqs, $seq;
      }
    @seqs = sort {length($b) <=> length($a)} @seqs;
    return @seqs;
  }

################################################ subroutine header begin ##

=head2 get_genomic_sequence_for_feature

 Usage     : my $seq = $genome_obj->get_genomic_sequence_for_feature($feature_obj)
 Purpose   : Gets the genomic sequence of a feature using a feature object, and of a user
           : defined sequence length before and after the feature
 Returns   : string (gemomic DNA of feature)
 Argument  : CoGe::Genome::DB::Feature object
 Throws    : none
 Comments  : The sequence is returned in 5'->3' orientation.
           : Also have check so that for features from datasets without genomic sequence
           : this will try to find an associated dataset with genomic sequence from which
           : to draw a genomic sequence for the feature.

See Also   : CoGe::Genome::DB::Data_information::get_associated_datasets

=cut

################################################## subroutine header end ##

sub get_genomic_sequence_for_feature
  {
    my $self = shift;
    my $feat;
    my ($upstream, $downstream) = (0,0);
    if (scalar @_ == 1) {$feat = shift;} #for legecy code
    else 
      {
	my %opts = @_;
	$feat = $opts{feature} || $opts{feat};
	$upstream = $opts{upstream} || $opts{up} || $opts{us} || 0;
	$downstream = $opts{downstream} || $opts{down} || $opts{ds} || 0;
      }
    my $seq;
    if ($feat->has_genomic_sequence())
      {
	foreach my $loc (sort {$a->start <=> $b->start} $feat->locs)
	  {
	    my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
								   start  => $loc->start,
								   stop   => $loc->stop,
								   chr    => $loc->chr,
								   org_id => $feat->data_info->organism->id,
								   info_id=> $feat->data_info->id,
								   strand => $loc->strand,
								  );
	    $seq .= $tmp_seq if $tmp_seq;
	  }

	if ($upstream)
	  {
	    my ($start, $stop);
	    if ($feat->strand =~ /-/)
	      {
		$stop = $feat->stop+$upstream;
		$start = $feat->stop+1;
	      }
	    else
	      {
		$start = $feat->start-$upstream;
		$stop = $feat->start-1;
	      }
	    my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
								   start  => $start,
								   stop   => $stop,
								   chr    => $feat->chr,
								   org_id => $feat->data_info->organism->id,
								   info_id=> $feat->data_info->id,
								   strand => $feat->strand,
								  );
	    $seq = $tmp_seq . $seq if $tmp_seq;

	  }
	if ($downstream)
	  {
	    my ($start, $stop);
	    if ($feat->strand =~ /-/)
	      {
		$start = $feat->start-$downstream;
		$stop = $feat->start-1;
	      }
	    else
	      {
		$stop = $feat->stop+$downstream;
		$start = $feat->stop+1;
	      }
	    my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
								   start  => $start,
								   stop   => $stop,
								   chr    => $feat->chr,
								   org_id => $feat->data_info->organism->id,
								   info_id=> $feat->data_info->id,
								   strand => $feat->strand,
								  );
	    $seq .= $tmp_seq if $tmp_seq;

	  }
      }
    else
      {
	outer: foreach my $ds ($feat->dataset->get_associated_datasets())
	  {
	    if ($ds->has_genomic_sequence)
	      {
		foreach my $chr ($ds->chromosomes)
		  {
		    if ($chr eq $feat->chr)
		      {
			foreach my $loc (sort {$a->start <=> $b->start} $feat->locs)
			  {
			    my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
										   start  => $loc->start,
										   stop   => $loc->stop,
										   chr    => $loc->chr,
										   org_id => $feat->data_info->organism->id,
										   info_id=> $ds->id,
										   strand => $loc->strand,
										  );
			    $seq .= $tmp_seq if $tmp_seq;
			  }
			if ($upstream)
			  {
			    my ($start, $stop);
			    if ($feat->strand =~ /-/)
			      {
				$stop = $feat->stop+$upstream;
				$start = $feat->stop+1;
			      }
			    else
			      {
				$start = $feat->start-$upstream;
				$stop = $feat->start-1;
			      }

			    my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
										   start  => $start,
										   stop   => $stop,
										   chr    => $feat->chr,
										   org_id => $feat->data_info->organism->id,
										   info_id=> $ds->id,
										   strand => $feat->strand,
										  );
			    $seq = $tmp_seq . $seq if $tmp_seq;
			    
			  }
			if ($downstream)
			  {
			    my ($start, $stop);
			    if ($feat->strand =~ /-/)
			      {
				$start = $feat->start-$downstream;
				$stop = $feat->start-1;
			      }
			    else
			      {
				$stop = $feat->stop+$downstream;
				$start = $feat->stop+1;
			      }
			    my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
										   start  => $start,
										   stop   => $stop,
										   chr    => $feat->chr,
										   org_id => $feat->data_info->organism->id,
										   info_id=> $ds->id,
										   strand => $feat->strand,
										  );
			    $seq .= $tmp_seq if $tmp_seq;
			    
			  }
  
			last outer;
		      }
		  }
	      }
	  }
      }
    return $seq;
  }

################################################ subroutine header begin ##

=head2 get_genomic_seq_for_feat

 Usage     : my $seq = $genome_obj->get_genomic_seq_for_feat($feature_obj)
 Purpose   : Alias for get_genomic_sequence_for_feature
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_genomic_seq_for_feat
  {
    my $self = shift;
    my $feat = shift;
    return $self->get_genomic_sequence_for_feature($feat);
  }


################################################ subroutine header begin ##

=head2 get_features_for_organism_and_version

 Usage     : my @feats = $db->get_features_for_organism_and_version(org=>$org,
                                                                    version=>6,
								    );
 Purpose   : Get all the features for an organism of a particular version
           : as specified by the data_information table.
 Returns   : an array or array_ref (wantarray) of CoGe::Genome::DB::Feature objects
 Argument  : org => either the name of the organism (which will be looked up in the database) --OR--
                    an CoGe::Genome:DB::Organism object
             version=> an integer that specifies the data inforation version from which to find features
 Throws    : Will return undef if no organism object can be found.

See Also   : 

=cut

################################################## subroutine header end ##


sub get_features_for_organism_and_version
  {
    my $self = shift;
    my %opts = @_;
    my $org = $opts{org} || $opts{organism};
    unless (ref ($org) =~ /organism/i)
      {
         ($org) = $self->get_org_obj->search_like(name=>"$org%");
      }
    return unless ref($org) =~ /organism/i;
    my $version = $opts{version} || $opts{ver} || $opts{v};
    my @feats;
    foreach my $di ($org->data_infos)
      {
         next unless $di->version == $version;
	 push @feats, $di->feats;
      }
    return wantarray ? @feats : \@feats;
  }

################################################ subroutine header begin ##

=head2 find_overlapping_features

 Usage     : my @overlap_feats = $db->find_overlapping_features(feat=>$feat);
 Purpose   : find overlapping features given a feature
 Returns   : an array or array ref of CoGe::Genome::DB::Feature objects
 Argument  : feat=> CoGe::Genome::DB::Feature object
             searchall => 1
                          searchs all datasets for the organism to 
                          which the feature belongs
 Throws    : 
 Comments  : If searchall is set, then all datasets for the organism to which the
             feature belongs are searched for overlapping features.  Otherwise, only 
             the dataset to which the feature belongs is used.

See Also   : CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##

sub find_overlapping_features
  {
    my $self = shift;
    my %opts = @_;
    my $feat = $opts{feat};
    return unless ref ($feat) =~ /Feat/i;
    my @di = $opts{searchall} ? $feat->dataset->organism->datasets : $feat->dataset;
    my @overlap_feats;
    foreach my $di (@di)
      {
	push @overlap_feats, $feat->get_features_in_region(start=>$feat->start,
							   stop =>$feat->stop,
							   chr  =>$feat->chr,
							   info_id=>$di,
							  );
      }
    return wantarray ? @overlap_feats : \@overlap_feats;
  }


################################################ subroutine header begin ##

=head2 get_features_for_organism

 Usage     : my (@feats) = $db->get_features_for_org(org=> "Arabidopsis",
                                                     type=>"gene",
                                                     chr=>1,
                                                     );
 Purpose   : get a list of feature objects for an organism based on user 
             defined criteria
 Returns   : an array or array ref of CoGe::Genome::DB::Feature objects
 Argument  : a hash of key->value pairs
             org           => CoGe::Genome::DB::Organism object or organism database 
                              id or organism name (required!)
             version | ver => which version of data to use (if not specified this
                            routine will find the most current version of the data
                            and use that)
             type        => the type of feature to return (e.g. "gene").  If this is
                            not specified it will return all features. (optional)
             chr | chromosome => chromosome name to limit search to just one 
                                 chromosome. (optional)
 Throws    : a warning will be thrown if the organism is invalid
 Comments  : this calls CoGe::Genome::DB::Data_information->get_features_for_organism

See Also   : CoGe::Genome::DB::Feature
             CoGe::Genome::DB::Data_information
=cut

################################################## subroutine header end ##


sub get_features_for_organism
  {
    my $self = shift;
    return $self->get_dataset_obj->get_features_for_organism(@_);
  }

################################################ subroutine header begin ##

=head2 all_feature_types

 Usage     : 
 Purpose   : returns all the feature type objects.  Useful for building search lists
 Returns   : an array of feature type objects
 Argument  : none
 Throws    : 
 Comments  : calls Feature_type->retrieve_all()

See Also   : CoGe::Genome::DB::Feature_type

=cut

################################################## subroutine header end ##

sub all_feature_types
  {
    my $self = shift;
    return $self->get_feature_type_obj->retrieve_all;
  }

################################################ subroutine header begin ##

=head2 all_orgs

 Usage     : 
 Purpose   : returns all the organism objects.  Useful for building search lists
 Returns   : an array of organism objects
 Argument  : none
 Throws    : 
 Comments  : calls Organism->retrieve_all()

See Also   : CoGe::Genome::DB::Organism

=cut

################################################## subroutine header end ##

sub all_orgs
  {
    my $self = shift;
    return $self->get_organism_obj->retrieve_all;
  }

################################################ subroutine header begin ##

=head2 get_organism

 Usage     : my $org = $coge->get_organism("Arabidopsis");
 Purpose   : gets an organism object for an organism by database id or name
 Returns   : CoGe::Genome::DB::Organism object
 Argument  : organism database id or name
 Throws    : none
 Comments  : calls Organism->resolve_organism

See Also   : CoGe::Genome::DB::Organism

=cut

################################################## subroutine header end ##

sub get_organism
  {
    my $self = shift;
    return $self->get_organism_obj->resolve_organism(@_);
  }

################################################ subroutine header begin ##

=head2 get_genomic_sequence

 Usage     : $object->get_sequence(start   => $start, 
                                   stop    => $stop, 
                                   chr     => $chr,
                                   dataset_id => $data_info->id());

 Purpose   : gets the genomic sequence for the specified conditions
 Returns   : a string (containing the genomic sequence)
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             dataset_id => data_information id in database (obtained from a
                        CoGe::Data_information object)
             strand  => 1 or -1.  Default 1.
                        if negative strand is requested, the complement
                        of the dna seq will be returned
 Throws    : undef if no sequence is obtained
 Comments  : Genomic_sequence->get_sequence

See Also   : CoGe::Genome::DB::Genomic_sequence

=cut

################################################## subroutine header end ##

sub get_genomic_sequence
  {
    my $self = shift;
    my %args = @_;
    my $gso =  CoGe::Genome::DB::Genomic_sequence->new();
    return $gso->get_sequence(%args);
  }

################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##



  

1; #this line is important and will help the module return a true value
