package CoGeX::Result::Experiment;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;
use base 'DBIx::Class::Core';
#use CoGeX::ResultSet::Experiment;
use File::Spec::Functions;
use Data::Dumper;
#use Text::Wrap;
use POSIX;
use Carp;
use LWP::Simple;

=head1 NAME

CoGeX::DatasetGroup

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<dataset_group> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<experiment_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<dataset_group_id> 
Type: INT, Default: undef, Nullable: no, Size: 11

C<data_source_id>
Type: INT, Default: 0, Nullable: 0, Size: 11

C<name>
Type:VARCHAR, Default: "", Nullable: no, Size: 255

C<description>
Type: Text, Default: undef, Nullable: yes

C<version>
Type:VARCHAR, Default: undef, Nullable: no, Size: 50

C<storage_path>
Type: VARCHAR, Default: undef, Nullable: 0, Size: 255

C<link>
Type: TEXT, Defaullt:  undef, Nullable: 1

C<restricted>
Type: INT, Default: 0, Nullable: 0, Size: 1




Belongs to CCoGeX::Result::DataSource> via C<data_source_id>
Belongs to CCoGeX::Result::DatasetGroup> via C<dataset_group_id>
Has many CCoGeX::Result::ExperimentTypeConnector> via C<experiment_id>
Has many CCoGeX::Result::ExperimentAnnotation> via C<experiment_id>


=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("experiment");
__PACKAGE__->add_columns(
  "experiment_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "dataset_group_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "data_source_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "description",
  {
    data_type => "TEXT",
    default_value => undef,
    is_nullable => 1,
  },
  "version",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 0,
    size => 50,
  },
  "storage_path",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 0,
    size => 255,
  },
  "restricted",  { data_type => "int", default_value => "0", is_nullable => 0, size => 1 },
  "access_count",  { data_type => "int", default_value => "0", is_nullable => 1, size => 10 },
  "link",
  {
    data_type => "text",
    default_value => undef,
    is_nullable => 1,
  },
);

__PACKAGE__->set_primary_key("experiment_id");
__PACKAGE__->has_many("experiment_type_connectors" => "CoGeX::Result::ExperimentTypeConnector", 'experiment_id');
__PACKAGE__->has_many("experiment_annotations" => "CoGeX::Result::ExperimentAnnotation", 'experiment_id');
__PACKAGE__->belongs_to("data_source" => "CoGeX::Result::DataSource", 'data_source_id');
__PACKAGE__->belongs_to("dataset_group" => "CoGeX::Result::DatasetGroup", 'dataset_group_id');

#Need something like UserGroups permission stuff in the future
#__PACKAGE__->has_many("user_group_data_connectors"=>"CoGeX::Result::UserGroupDataConnector","dataset_group_id");

################################################ subroutine header begin ##

=head2 desc

 Usage     : 
 Purpose   : alias for $self->description
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub desc
  {
    return shift->description(@_);
  }


################################################ subroutine header begin ##

=head2 source

 Usage     : 
 Purpose   : alias for $self->data_source
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub source
  {
    shift->data_source(@_);
  }

################################################ subroutine header begin ##

=head2 genome

 Usage     : 
 Purpose   : alias for $self->dataset_group
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub genome
  {
    shift->dataset_group(@_);
  }

################################################ subroutine header begin ##

=head2 annotations

 Usage     : 
 Purpose   : alias for $self->experiment_annotations
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub annotation
  {
    shift->experiment_annotations(@_);
  }


################################################ subroutine header begin ##

=head2 experiment_types

 Usage     : $self->experiment_types
 Purpose   : pass through experiment_type_connector to fake a many-to-many connection with experiment_type
 Returns   : 
 Argument  : none, really
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub experiment_types
 {
   map {$_->experiment_type} shift->experiment_type_connectors();
 }


################################################ subroutine header begin ##

=head2 type

 Usage     : 
 Purpose   : alias for $self->experiment_types
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub type
  {
    shift->experiment_type(@_);
  }



=head1 BUGS

Whole lotta

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

=cut

1;
