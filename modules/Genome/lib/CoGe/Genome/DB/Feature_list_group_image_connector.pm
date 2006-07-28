package CoGe::Genome::DB::Feature_list_group_image_connector;
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
    __PACKAGE__->set_up_table('feature_list_group_image_connector');
    __PACKAGE__->has_a(feature_list_group_id=>'CoGe::Genome::DB::Feature_list_group');
    __PACKAGE__->has_a(image_id=>'CoGe::Genome::DB::Image');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature_list_group_image_connector

=head1 SYNOPSIS

=head1 DESCRIPTION

This object connects the feature_list_group table eithe image table.  



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

 feature_list_group_id => returns the associate Coge::Genome::DB::Feature_list_group 
                          object 
 feature_list_group
 flg
 group

 image_id => returns the associated CoGe::Genome::DB::Image object
 image

 feature_list_group_image_connector_id =>database id
 id

 

 new              =>  creates a new object (inherited from Class::Accessor)

=cut


sub flg
  {
    my $self = shift;
    return $self->feature_list_group_id;
  }

sub group
  {
    my $self = shift;
    return $self->feature_list_group_id;
  }

sub image
  {
    my $self = shift;
    return $self->image_id;
  }

sub id
  {
    my $self = shift;
    return $self->feature_list_group_image_connector_id();
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

