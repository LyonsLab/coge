package CoGe::Genome::DB::Feature_list;
use strict;
use base 'CoGe::Genome::DB';

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->set_up_table('feature_list');
    __PACKAGE__->has_many(feature_list_connectors=>'CoGe::Genome::DB::Feature_list_connector');
    __PACKAGE__->has_many(user_group_feature_list_permission_connectors=>'CoGe::Genome::DB::User_group_feature_list_permission_connector');
    __PACKAGE__->has_a('feature_list_group_id' => 'CoGe::Genome::DB::Feature_list_group');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature_list

=head1 SYNOPSIS

=head1 DESCRIPTION

This object accesses the feature_list table.  The feature_list contains the name and
description of a feature_list.  Connections between feature_lists and features occurs
in the feature_list_connector table.

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

CoGe::Genome
CoGe::Genome::DB
Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

 name             =>  name of feature list

 description      =>  description
 desc             =>  alias for description

 feature_list_id  =>  database entry id
 id               =>  alias for feature_list_id

 feature_list_group_id => returns the CoGe::Genome::DB::Feature_list_group object
                          associated with the feature_list (if one exists)
 feature_list_group
 fl_group
 group
 flg  

 feature_list_connectors
 feature_list_connector
 fl_connector
 flc

 user_group_feature_list_permission_connectors => 
 user_group_feature_list_permission_connector  => 
 ugflp_connector
 ugflpc

 new              =>  creates a new object (inherited from Class::Accessor)

=cut


sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub feature_list_group
  {
    my $self = shift;
    return $self->feature_list_group_id();
  }

sub fl_group
  {
    my $self = shift;
    return $self->feature_list_group_id();
  }

sub group
  {
    my $self = shift;
    return $self->feature_list_group_id();
  }

sub flg
  {
    my $self = shift;
    return $self->feature_list_group_id();
  }

sub id
  {
    my $self = shift;
    return $self->feature_list_id();
  }

sub feature_list_connector
  {
    my $self = shift;
    return $self->feature_list_connectors();
  }

sub fl_connector
  {
    my $self = shift;
    return $self->feature_list_connectors();
  }

sub flc
  {
    my $self = shift;
    return $self->feature_list_connectors();
  }

sub user_group_feature_list_permission_connector
  {
    my $self = shift;
    return $self->user_group_feature_list_permission_connectors();
  }

sub ugflp_connector
  {
    my $self = shift;
    return $self->user_group_feature_list_permission_connectors();
  }

sub ugflpc
  {
    my $self = shift;
    return $self->user_group_feature_list_permission_connectors();
  }

################################################ subroutine header begin ##

=head2 features

 Usage     : my @features = $feature_list_obj->features();
 Purpose   : fetches the feature objects assoicated with a feature list
             by hopping though the feature_list_connector table
 Returns   : an array or array ref depending on wantarray
 Argument  : none
 Throws    : none
 Comments  : 

See Also   : CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##

sub features
  {
    my $self = shift;
    my @features;
    foreach my $flc ($self->flc)
      {
	push @features, $flc->feature;
      }
    return wantarray ? @features : \@features
  }

sub readers
  {
    my $self = shift;
    my @readers;
    foreach my $ugflpc ($self->ugflpc)
      {
	foreach my $perm ($ugflpc->permissions)
	  {
	    if ($perm->name =~ /feature list read/)
	      {
		push @readers, $ugflpc->user_group;
	      }
	  }
      }
    return wantarray ? @readers : \@readers;
  }

sub editors
  {
    my $self = shift;
    my @editors;
    foreach my $ugflpc ($self->ugflpc)
      {
	foreach my $perm ($ugflpc->permissions)
	  {
	    if ($perm->name =~ /feature list edit/)
	      {
		push @editors, $ugflpc->user_group;
	      }
	  }
      }
    return wantarray ? @editors : \@editors;
  }

################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

1; #this line is important and will help the module return a true value

