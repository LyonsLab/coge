package CoGe::Genome::DB::Permission;
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
    __PACKAGE__->set_up_table('permission');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Permission

=head1 SYNOPSIS

=head1 DESCRIPTION

 This object accesses the information in the "permission" table of the CoGe::Genome database.
 This table has two columns -- name and description.  This structure allows for the specification 
 and expansion of any permission type that can then be used to check to see what a user can do (see user
 and user_group tables and objects).  Examples of permissions are:

 Name                    Description
 admin                   Full administrative privileges
 read feature list       Read a feature list
 edit feature list       Add/Delete features from a feature list
 create feature list     Create and delete a feature list
 

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

 name             =>  name of permission

 description      =>  description
 desc             =>  alias for description

 permission_id    =>  database entry id
 id               =>  alias for permission_id

 

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
    return $self->permission_id();
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

