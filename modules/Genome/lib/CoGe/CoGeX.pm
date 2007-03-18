package CoGe::CoGeX;


use strict;
use Data::Dumper;
use warnings;
use base 'DBIx::Class::Schema';

__PACKAGE__->load_classes qw(
Feature
Annotation
FeatureList
Image
UserGroup
AnnotationType
FeatureListConnector
Location
UserGroupConnector
AnnotationTypeGroup
FeatureListGroup
Organism
UserGroupFeatureListPermissionConnector
FeatureListGroupImageConnector
Permission
UserSession
DataSource
FeatureName
Sequence
Dataset
FeatureType
SequenceType
GenomicSequence
User
);


#use base 'CoGeX';

use CoGe::CoGeX::Feature;
use CoGe::CoGeX::Annotation;
use CoGe::CoGeX::FeatureList;
use CoGe::CoGeX::Image;
use CoGe::CoGeX::UserGroup;
use CoGe::CoGeX::AnnotationType;
use CoGe::CoGeX::FeatureListConnector;
use CoGe::CoGeX::Location;
use CoGe::CoGeX::UserGroupConnector;
use CoGe::CoGeX::AnnotationTypeGroup;
use CoGe::CoGeX::FeatureListGroup;
use CoGe::CoGeX::Organism;
use CoGe::CoGeX::UserGroupFeatureListPermissionConnector;
use CoGe::CoGeX::FeatureListGroupImageConnector;
use CoGe::CoGeX::Permission;
use CoGe::CoGeX::UserSession;
use CoGe::CoGeX::DataSource;
use CoGe::CoGeX::FeatureName;
use CoGe::CoGeX::Sequence;
use CoGe::CoGeX::Dataset;
use CoGe::CoGeX::FeatureType;
use CoGe::CoGeX::SequenceType;

use CoGe::CoGeX::GenomicSequence;
use CoGe::CoGeX::User;

use vars qw( $VERSION );

$VERSION = 0.01;

=head1 NAME

CoGeX - CoGeX

=head1 SYNOPSIS

  use CoGeX;
  blah blah blah


=head1 DESCRIPTION

Primary object for interacting with CoGe database system.

=head1 USAGE

  use CoGeX;

  my $connstr = 'dbi:mysql:genomes:biocon:3306';
  my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' ); # Biocon's ro user

  my $rs = $s->resultset('Feature')->search(
                {
                  'organism.name' => "Arabidopsis thaliana"
                },
                { join => [ 'dataset', 'organism' ] }
  );

=head1 BUGS


=head1 SUPPORT


=head1 AUTHOR

    Brian C. Thomas
    CPAN ID: BCT
    UC - Berkeley
    bcthomas@nature.berkeley.edu
    http://toxic.berkeley.edu

    Brent Pedersen
    CPAN ID: BPEDERSE
    UC - Berkeley
    bpederse@nature.berkeley.edu
    http://toxic.berkeley.edu


    Eric Lyons
    CPAN ID:
    UC - Berkeley
    elyons@nature.berkeley.edu
    http://toxic.berkeley.edu

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

  sub annotation
    {
      return new CoGe::CoGeX::Annotation;
    }

  sub featureList
    {
      return new CoGe::CoGeX::FeatureList;
    }

  sub image
    {
      return new CoGe::CoGeX::Image;
    }

  sub userGroup
    {
      return new CoGe::CoGeX::UserGroup;
    }
  sub annotationType
    {
      return new CoGe::CoGeX::AnnotationType;
    }
  sub featureListConnector
    {
      return new CoGe::CoGeX::FeatureListConnector;
    }
  sub location
    {
      return new CoGe::CoGeX::Location;
    }
  sub userGroupConnector
    {
      return new CoGe::CoGeX::UserGroupConnector;
    }
  sub anntationTypeGroup
    {
      return new CoGe::CoGeX::AnnotationTypeGroup;
    }
  sub featureListGroup
    {
      return new CoGe::CoGeX::FeatureListGroup;
    }
  sub organisn
    {
      return new CoGe::CoGeX::Organism;
    }
  sub userGroupFeatureListPermisionConnector
    {
      return new CoGe::CoGeX::UserGroupFeatureListPermissionConnector;
    }
  sub featureListGroupImageConnector
    {
      return new CoGe::CoGeX::FeatureListGroupImageConnector;
    }
  sub permission
    {
      return new CoGe::CoGeX::Permission;
    }
  sub userSession
    {
      return new CoGe::CoGeX::UserSession;
    }
  sub dataSource
    {
      return new CoGe::CoGeX::DataSource;
    }
  sub featureName
    {
      return new CoGe::CoGeX::FeatureName;
    }
  sub sequence
    {
      return new CoGe::CoGeX::Sequence;
    }
  sub dataset
    {
      return new CoGe::CoGeX::Dataset;
    }
  sub featureType
    {
      return new CoGe::CoGeX::FeatureType;
    }
  sub sequenceType
    {
      return new CoGe::CoGeX::SequenceType;
    }
  sub feature
    {
      return new CoGe::CoGeX::Feature;
    }
  sub genomicSequence
    {
      return new CoGe::CoGeX::GenomicSequence;
    }
  sub user
    {
      return new CoGe::CoGeX::User;
    }

 sub reverse_complement
   {
     my $self = shift;
     print STDERR Dumper \@_;
     $self->feature->reverse_complement(@_);
   }

################################################ subroutine header begin ##

=head2 get_features_in_region

 Usage     : $object->get_features_in_region(start   => $start, 
                                             stop    => $stop, 
                                             chr     => $chr,
                                             dataset_id => $dataset->id());

 Purpose   : gets all the features in a specified genomic region
 Returns   : an array or an array_ref of feature objects (wantarray)
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             dataset_id => dataset id in database (obtained from a
                        CoGe::Dataset object)
                        of the dna seq will be returned
             OPTIONAL
             count   => flag to return only the number of features in a region
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_features_in_region
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{'start'} || $opts{'START'} || $opts{begin} || $opts{BEGIN};
    $start = 0 unless $start;
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    $stop = $start unless defined $stop;
    my $chr = $opts{chr} || $opts{CHR} || $opts{chromosome} || $opts{CHROMOSOME};
    my $dataset_id = $opts{dataset} || $opts{dataset_id} || $opts{info_id} || $opts{INFO_ID} || $opts{data_info_id} || $opts{DATA_INFO_ID} ;
    my $count_flag = $opts{count} || $opts{COUNT};
    if ($count_flag)
      {
	return $self->resultset('Feature')->count({
						   "locations.chromosome"=>$chr,
						   "dataset_id"=>$dataset_id,
						   "locations.stop"=>  {"<=" => $stop},
						   "locations.start"=> {">=" => $start},
						  },
						  {
						   join => ["locations"],
						   distinct=>['feature_id'],
						  }
						 );
      }
    my @feats = $self->resultset('Feature')->search({
						     "locations.chromosome"=>$chr,
						     "dataset_id"=>$dataset_id,
						     "locations.stop"=>  {"<=" => $stop},
						     "locations.start"=> {">=" => $start},
						    },
						    {
						     join => ["locations"],
						     distinct=>['feature_id'],						    }
						   );
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 count_features_in_region

 Usage     : $object->count_features_in_region(start   => $start, 
                                             stop    => $stop, 
                                             chr     => $chr,
                                             dataset_id => $dataset->id());

 Purpose   : counts the features in a specified genomic region
 Returns   : an integer
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             dataset_id => dataset id in database (obtained from a
                        CoGe::Dataset object)
                        of the dna seq will be returned
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub count_features_in_region
  {
    my $self = shift;
    my %opts = @_;
    return $self->get_features_in_region (%opts, count=>1);
  }


1;
