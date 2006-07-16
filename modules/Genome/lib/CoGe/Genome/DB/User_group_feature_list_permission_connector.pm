package CoGe::Genome::DB::User_group_feature_list_permission_connector;
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
    __PACKAGE__->set_up_table('user_group_feature_list_permission_connector');
    __PACKAGE__->has_a('user_group_id'=>'CoGe::Genome::DB::User_group');
    __PACKAGE__->has_a('feature_list_id'=>'CoGe::Genome::DB::Feature_list');
    __PACKAGE__->has_a('permission_id'=>'CoGe::Genome::DB::Permission');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::User_group_feature_list_permission_connector

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
Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

 user_group_id            => user_group table database id
 ug_id
 ugid
 user_group
 group
 ug

 feature_list_id          => user_group table database id
 fl_id
 flid
 feature_list
 feat_list
 flist
 list
 fl

 permission_id            => permission table database id
 p_id
 pid
 permissions
 permission
 perms
 perm
 p

 user_group_feature_list_permission_connector_id =>  database entry id
 id                       =>  alias for user_group_permission_id

 new                      =>  creates a new object (inherited from Class::Accessor)

=cut


sub ug_id
  {
    my $self = shift;
    return $self->user_group_id(@_);
  }

sub ugid
  {
    my $self = shift;
    return $self->user_group_id(@_);
  }

sub user_group
  {
    my $self = shift;
    return $self->user_group_id(@_);
  }

sub group
  {
    my $self = shift;
    return $self->user_group_id(@_);
  }

sub ug
  {
    my $self = shift;
    return $self->user_group_id(@_);
  }

sub fl_id
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }

sub flid
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }

sub feature_list
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }

sub feat_list
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }

sub flist
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }

sub list
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }

sub fl
  {
    my $self = shift;
    return $self->feature_list_id(@_);
  }



sub p_id
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub pid
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub permissions
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub permission
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub perms
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub perm
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub p
  {
    my $self = shift;
    return $self->permission_id(@_);
  }

sub id
  {
    my $self = shift;
    return $self->user_group_feature_list_permission_connector_id();
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

