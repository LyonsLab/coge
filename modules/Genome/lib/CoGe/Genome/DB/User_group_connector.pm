package CoGe::Genome::DB::User_group_connector;
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
    __PACKAGE__->set_up_table('user_group_connector');
    __PACKAGE__->has_a('user_id'=>'CoGe::Genome::DB::User');
    __PACKAGE__->has_a('user_group_id'=>'CoGe::Genome::DB::User_group');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::User_group_connector

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

 user_id                 => user table database id
 suer
 u_id          
 uid
 u

 user_group_id           => user group table database id
 user_group
 ug_id
 ugid
 ug

 user_group_connector_id =>  database entry id
 id                      =>  alias for user_group_connector_id

 new                     =>  creates a new object (inherited from Class::Accessor)

=cut


sub u_id
  {
    my $self = shift;
    return $self->user_id(@_);
  }

sub user
  {
    my $self = shift;
    return $self->user_id(@_);
  }

sub u
  {
    my $self = shift;
    return $self->user_id(@_);
  }

sub uid
  {
    my $self = shift;
    return $self->user_id(@_);
  }

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

sub ug
  {
    my $self = shift;
    return $self->user_group_id(@_);
  }

sub id
  {
    my $self = shift;
    return $self->user_group_connector_id();
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

