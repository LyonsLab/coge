package CoGe::Genome::DB::Feature_list_connector;
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
    __PACKAGE__->set_up_table('feature_list_connector');
    __PACKAGE__->has_a('feature_id'=>'CoGe::Genome::DB::Feature');
    __PACKAGE__->has_a('feature_list_id'=>'CoGe::Genome::DB::Feature_list');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature_list_connector

=head1 SYNOPSIS

=head1 DESCRIPTION

This object accesses the table which connects the feature table to the feature_list table.
This allows multiple features to be in multiple feature lists.

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

 user_group_id    =>  database entry id for the user_group table
 ug_id
 ugid               

 feature_list_id  =>  get the feature list object
 fl_id
 flid
 feature_id       =>  get the feature object
 feature
 feat

 permission_id    =>  database entry id for the permission table
 p_id
 pid

 preferred_name   =>  optional name that is preferred for the feature
 name             

 description          => user defined description of the feature
 desc


 feature_list_connector_id  =>  database entry id
 id                         =>  alias for feature_list_id

 new              =>  creates a new object (inherited from Class::Accessor)

=cut


sub ug_id
  {
    my $self = shift;
    return $self->user_group_id();
  }

sub ugid
  {
    my $self = shift;
    return $self->user_group_id();
  }

sub fl_id
  {
    my $self = shift;
    return $self->feature_list_id();
  }

sub feature
  {
    my $self = shift;
    return $self->feature_id();
  }

sub feat
  {
    my $self = shift;
    return $self->feature_id();
  }

sub p_id
  {
    my $self = shift;
    return $self->permission_id();
  }

sub pid
  {
    my $self = shift;
    return $self->permission_id();
  }

sub name
  {
    my $self = shift;
    return $self->preferred_name();
  }

sub desc
  {
    my $self = shift;
    return $self->description();
  }

sub id
  {
    my $self = shift;
    return $self->feature_list_connector_id();
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

