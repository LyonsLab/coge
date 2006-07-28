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
    __PACKAGE__->set_sql('fetch_ordered_lists' => qq{
SELECT *
  FROM feature_list fl
 ORDER BY ? ?
 LIMIT ?, ?
});
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

=head2 get_lists

 Usage     : my @lists = $obj->get_lists(start=>1, limit=>50, sort=>"name", sort_type=>"desc");
 Purpose   : Get feature lists from the database
 Returns   : returns an array or arrayref of CoGe::Genome::DB::Feature_list objects
 Argument  : hash of key-value pairs:
     start=> where in the list to start (default 1)
     limit=> number of items to retrieve from the database (default 200)
     sort => which column in the database on which to sort (default feature_list_id)
sort_type => ASC or DESC for ascending or descending sorts

 Throws    : 
 Comments  : If you want all the lists in the database use $obj->retrieve_all()

See Also   : 

=cut

################################################## subroutine header end ##



sub get_lists
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin};
    my $num = $opts{limit} || $opts{limit};
    my $sort = $opts{sort} || 'feature_list_id';
    my $sort_type = $opts{sort_type} || 'desc';
    $start = 1 unless $start;
    $start = $start -1 if $start;
    $num = 200 unless $num;
#    my $user = $opts{user}; will need this at some point
#    my $user_group = $opts{user_group};
    my $sth = $self->sql_fetch_ordered_lists;
    $sth->execute($sort, $sort_type, $start, $num);
    my @fls;
    while (my $q = $sth->fetchrow_arrayref)
      {
	push @fls, $self->retrieve($q->[0]);
      }
    $sth->finish;
    return wantarray ? @fls : \@fls;
  }
################################################ subroutine header begin ##

=head2 get_preferred_names

 Usage     : my @names = $feature_list_obj->get_preferred_names();
 Purpose   : gets the perferred name for the feature list if they were 
             specified in the database
 Returns   : an array or array ref depending on wantarray
 Argument  : none
 Throws    : none
 Comments  : This information is stored in the feature_list_connector table

See Also   : CoGe::Genome::DB::Feature_list_connector

=cut

################################################## subroutine header end ##

sub get_preferred_names
  {
    my $self = shift;
    my @names ;
    foreach my $flc ($self->flc)
      {
	push @names, $flc->preferred_name if $flc->preferred_name;
      }
    return wantarray ? @names : \@names
 
  }

################################################ subroutine header begin ##

=head2 preferred_name

 Usage     : my $name = $feature_list_obj->preferred_name($feat_obj);
 Purpose   : find the preferred name for a feature as stored in the feature list connector
 Returns   : a string or undef
 Argument  : a feature_obj of a feature from the feature list
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub preferred_name
  {
    my $self = shift;
    my $feat = shift;
    foreach my $flc ($self->flc)
      {
	return $flc->name if $flc->feature_id eq $feat->id;
      }
    return;
  }
################################################ subroutine header begin ##

=head2 preferred_description

 Usage     : my $desc = $feature_list_obj->preferred_description($feat_obj);
 Purpose   : find the description for the preferred name for a feature as stored in the feature list connector
 Returns   : a string or undef
 Argument  : a feature_obj of a feature from the feature list
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub preferred_description
  {
    my $self = shift;
    my $feat = shift;
    foreach my $flc ($self->flc)
      {
	return $flc->desc if $flc->feature_id eq $feat->id;
      }
    return;
  }

sub preferred_desc
  {
    my $self = shift;
    return $self->preferred_description(@_);
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

