package CoGe::Genome::DB::Feature_list_group;
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
    __PACKAGE__->set_up_table('feature_list_group');
    __PACKAGE__->has_many(feature_lists=>'CoGe::Genome::DB::Feature_list');
    __PACKAGE__->has_many(feature_list_group_image_connectors=>'CoGe::Genome::DB::Feature_list_group_image_connector');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature_list_group

=head1 SYNOPSIS

=head1 DESCRIPTION

This object accesses the feature_list_group table.  This is to organize single
or multiple feature_lists into a larger group.


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

 name             =>  name of feature list group

 description      =>  description
 desc             =>  alias for description

 feature_list_group_id  =>  database entry id
 id               =>  alias for feature_list_id

 feature_lists    => returns Coge::Genome::DB::Feature_list objects associated
                     with the feature_list_group
 feature_list
 lists
 list
 fl

 feature_list_group_image_connectors => returns Feature_list_group_image_conector
                                        objects
 flgi_conenctor
 flgic

 images => returns Image objects linked through the feature_list_group_image_connector
           object
 image

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
    return $self->feature_list_group_id();
  }

sub feature_list
  {
    my $self = shift;
    return $self->feature_lists();
  }

sub lists
  {
    my $self = shift;
    return $self->feature_lists();
  }

sub list
  {
    my $self = shift;
    return $self->feature_lists();
  }

sub fl
  {
    my $self = shift;
    return $self->feature_lists();
  }

sub flgi_connector
  {
    my $self = shift;
    return $self->feature_list_group_image_connectors();
  }

sub flgic
  {
    my $self = shift;
    return $self->feature_list_group_image_connectors();
  }

sub images
  {
    my $self = shift;
    my @images = map{$_->image} $self->flgic;
  }

sub image
  {
    my $self = shift;
    return $self->images;
  }

sub insert_image
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $img = $opts{img};
    my $db = new CoGe::Genome;
    my $io = $db->get_image_obj->insert({
					 name=>$name,
					 description=>$desc,
					 image=>$img,
					});
    my $ico = $db->get_feature_list_group_image_connector_obj->insert({
								   feature_list_group_id=>$self->id, 
								   image_id=>$io->id,
								  });
  }

################################################ subroutine header begin ##

=head2 features

 Usage     : my @features = $feature_list_group_obj->features();
 Purpose   : fetches the feature objects assoicated with a feature list group
             by hopping though the feature_list_connector table
 Returns   : an array or array ref depending on wantarray
 Argument  : none
 Throws    : none
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_list
             CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##



sub features
  {
    my $self = shift;
    my %feats;
    foreach my $fl ($self->fl)
      {
	foreach my $f ($fl->features)
	  {
	    $feats{$f->id} = $f;
	  }
      }
    return wantarray ? values %feats : [values %feats];
  }

################################################ subroutine header begin ##

=head2 preferred_names

 Usage     : my @names = $feature_list_group obj->preferred_names();
 Purpose   : gets the perferred name for the feature list group if they were 
             specified in the database
 Returns   : an array or array ref depending on wantarray
 Argument  : none
 Throws    : none
 Comments  : This information is stored in the feature_list_connector table

See Also   : CoGe::Genome::DB::Feature_list_connector
             CoGe::Genome::DB::Feature_list

=cut

################################################## subroutine header end ##

sub preferred_names
  {
    my $self = shift;
    my %seen;
    my @names = grep { !$seen{$_} ++} map {$_->get_preferred_names} $self->fl;
    return wantarray ? @names : \@names;
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

