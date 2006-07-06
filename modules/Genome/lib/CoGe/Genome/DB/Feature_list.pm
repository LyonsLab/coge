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
    return $self->feature_list_id();
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

