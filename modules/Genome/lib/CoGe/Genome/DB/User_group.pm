package CoGe::Genome::DB::User_group;
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
    __PACKAGE__->set_up_table('user_group');
    __PACKAGE__->has_many(user_group_connectors=>'CoGe::Genome::DB::User_group_connector');
    __PACKAGE__->has_many(user_group_feature_list_permission_connectors=>'CoGe::Genome::DB::User_group_feature_list_permission_connector');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::User_group

=head1 SYNOPSIS

=head1 DESCRIPTION

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
CoGe::Genome::DB::Sequence
Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

 name             =>  name of user group

 description      =>  description
 desc             =>  alias for description

 user_group_id          =>  database entry id
 id               =>  alias for sequence_type_id

 user_group_conenctors
 user_group_connector
 ug_connector
 ugc

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

sub id
  {
    my $self = shift;
    return $self->user_group_id();
  }

sub user_group_connector
  {
    my $self = shift;
    return $self->user_group_connectors();
  }

sub ug_connector
  {
    my $self = shift;
    return $self->user_group_connectors();
  }

sub ugc
  {
    my $self = shift;
    return $self->user_group_connectors();
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

