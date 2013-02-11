package CoGeX::Result::User;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';
use Data::Dumper;

=head1 NAME

CoGeX::User

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<user_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<user_name>
Type: VARCHAR, Default: "", Nullable: no, Size: 10

C<first_name>
Type: VARCHAR, Default: "", Nullable: no, Size: 10

C<last_name>
Type: VARCHAR, Default: "", Nullable: no, Size: 10

C<email>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 50

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255


Has many CCoGeX::Result::UserSession> via C<user_id>

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("user");
__PACKAGE__->add_columns(
	"user_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
	"user_name",
	{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 10 },
	"first_name",
	{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 10 },
	"last_name",
	{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 10 },
	"email",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 50 },
	"description",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 255 },
	"image_id",
	{ data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
	"date",
	{ data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
	"admin",
	{ data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 1 }
);
__PACKAGE__->set_primary_key("user_id");
__PACKAGE__->has_many( 'sessions'  => "CoGeX::Result::UserSession", 'user_id' );
__PACKAGE__->has_many( 'works'     => "CoGeX::Result::Work",        'user_id' );
__PACKAGE__->has_many( 'workflows' => "CoGeX::Result::Workflow",    'user_id' );
#__PACKAGE__->has_many( 'user_group_connectors' => "CoGeX::Result::UserGroupConnector", 'user_id' );
__PACKAGE__->has_many( 'user_connectors' => "CoGeX::Result::UserConnector", { "foreign.parent_id" => "self.user_id" } );
__PACKAGE__->has_many( 'logs' => "CoGeX::Result::Log",    'user_id' );
__PACKAGE__->has_many( 'jobs' => "CoGeX::Result::Job",    'user_id' );
__PACKAGE__->belongs_to( image => 'CoGeX::Result::Image', 'image_id');


################################################ subroutine header begin ##

=head2 generate_passwd

 Usage     : 
 Purpose   : Generates a password based on a hashed string and a salt value.
 Returns   : Hash of password and salt value.
 Argument  : 'passwd' or 'pwd'
 Throws    : None
 Comments  : 

See Also   : check_passwd()

=cut

################################################## subroutine header end ##

sub generate_passwd {
	my $self      = shift;
	my %opts      = @_;
	my $pwd       = $opts{passwd} || $opts{pwd};
	my $crypt_pwd = crypt( $pwd, "12" );
}

################################################ subroutine header begin ##

=head2 check_passwd

 Usage     : 
 Purpose   : Checks to see if entered password matches user password.
 Returns   : Result of logic test 'eq' between password has from the database and a hash of the user supplied password.
 Argument  : 'passwd' or 'pwd'
 Throws    : None
 Comments  : 

See Also   : generate_passwd()

=cut

################################################## subroutine header end ##

sub check_passwd {
	my $self = shift;
	my %opts = @_;
	my $pwd  = $opts{passwd} || $opts{pwd};
	return crypt( $pwd, $self->passwd ) eq $self->passwd;
}

################################################ subroutine header begin ##

=head2 name

 Usage     : 
 Purpose   : alias for $self->user_name
 Returns   : string
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub name {
	my $self = shift;
	return $self->user_name(@_);
}

sub display_name {
	my $self = shift;
	my $name = $self->user_name;
	$name = $self->first_name if $self->first_name;
	$name .= ' ' . $self->last_name if $self->first_name && $self->last_name;
	return $name;	
}

################################################ subroutine header begin ##

=head2 user_groups

 Usage     : 
 Purpose   : Returns the set of groups a user belongs to
 Returns   : Array of Groups
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub user_groups {
	my $self = shift;
	return unless $self->id; # ignore public user

	# my @user_groups = ();
	# foreach my $user_group_connector ( $self->user_group_connectors() ) {
	# 	push( @user_groups, $user_group_connector->user_group() )
	# 	  if $user_group_connector->user_group();
	# }

	my @groups;
	foreach my $conn ($self->user_connectors( {parent_type=>$node_types->{user}, child_type=>$node_types->{group} }) ) {
		push @groups, $conn->child_group;
	}
	return wantarray ? @groups : \@groups;
}

sub groups {
	return shift->user_groups(@_);
}

################################################ subroutine header begin ##

=head2 collaborators

 Usage     : 
 Purpose   : return user's collaborators
 Returns   : wantarray of user objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub collaborators {
	my $self = shift;

	my %users;
	foreach my $group ( $self->groups ) {
		foreach my $user ( $group->users ) {
			next if ($user->id == $self->id);
			$users{$user->id} = $user;
		}
	}

	return wantarray ? values %users : [ values %users ];
}

sub has_collaborator {
	my $self = shift;
	my $user = shift;
	my $uid = $user =~ /^\d+$/ ? $user : $user->id;
	foreach my $group ( $self->groups ) {
		foreach my $user ( $group->users ) {
			return 1 if ($user->id == $uid);
		}
	}	
	return 0;
}

################################################ subroutine header begin ##

=head2 owner_group

 Usage     : 
 Purpose   : return user's owner group
 Returns   : user_group object
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

# sub owner_group {
# 	my $self = shift;
# 	foreach my $group ( $self->groups ) {
# 		return $group if ($group->locked && $group->role->name =~ /owner/i);
# 	}
# 	return;
# }

################################################ subroutine header begin ##

=head2 owner_list

 Usage     : return user's owner list
 Purpose   : 
 Returns   : list object
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

# sub owner_list {
# 	return shift->owner_group->owner_list;
# }

################################################ subroutine header begin ##

=head2 shared_group

 Usage     : 
 Purpose   : return user's shared group
 Returns   : user_group object
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

# sub shared_group { #FIXME needs to be removed with addition of user_connector table
# 	my $self = shift;
# 	foreach my $group ( $self->groups ) {
# 		return $group if ($group->locked && $group->name eq $self->name && $group->description =~ /shared/i && $group->role->name =~ /reader/i);
# 	}
# 	return;
# }

################################################ subroutine header begin ##

=head2 owner_list

 Usage     : return user's shared list
 Purpose   : 
 Returns   : list object
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

# sub shared_list { #FIXME needs to be removed with addition of user_connector table
# 	my $self = shift;
# 	my $group = $self->shared_group;
# 	return if (not $group);
# 	return $group->shared_list;
# }

################################################ subroutine header begin ##

=head2 is_admin

 Usage     : $self->is_admin
 Purpose   : determine if a user is an admin
 Returns   : 1 if an admin, 0 if not
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub is_admin {
	# my $self = shift;
	# foreach my $group ( $self->groups ) {
	# 	return 1 if $group->role->name =~ /admin/i;
	# }
	# return 0;
	return shift->admin;
}

################################################ subroutine header begin ##

=head2 has_access_to_...

 Usage     : 
 Purpose   : checks to see if a user has access to a ...
 Returns   : 1/0
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub has_access_to_list {
	my $self = shift @_;
	my $list  = shift;
	return $self->has_access(list => $list);
}

sub has_access_to_genome {
	my $self = shift @_;
	my $dsg  = shift;
	return $self->has_access(dsg => $dsg);
}

sub has_access_to_experiment {
	my $self = shift @_;
	my $experiment = shift;
	return $self->has_access(experiment => $experiment);
}

sub has_access_to_dataset {
	my $self = shift @_;
	my $ds   = shift;
	return $self->has_access(ds => $ds);
}

################################################ subroutine header begin ##

=head2 has_access

 Usage     : 
 Purpose   : checks to see if a user has access to a dataset/genome/list/experiment
 Returns   : 1/0
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub has_access {
	my $self = shift;
	my %opts = @_;
	my $dsg  = $opts{dsg};
	my $ds   = $opts{ds};
	my $list = $opts{list};
	my $experiment = $opts{experiment};
	return 0 unless $dsg || $ds || $list || $experiment;
	return 0 unless $self->id;

	# foreach my $group ( $self->groups ) {
	# 	return 1 if $group->role->name eq "Admin";
	# }
	return 1 if $self->is_admin;

	if ($ds) {
		my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
		foreach my $ds ( $self->datasets() ) {
			if ( $dsid == $ds->id ) {
				return 1;
			}
		}
	}
	
	if ($dsg) {
		my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
		my %genomes;

		#get gids directly connected to user
#		map {$genomes{$_->child_id}=1} $self->user_connectors({child_type=>2});
		#check user lists to genomes
#		foreach my $list ($self->lists)
#		  {
		    
#		  }

		foreach my $genome ( $self->genomes(include_deleted => 1, ids=>1) ) {
			if ( $dsgid == $genome ) {
				return 1;
			}
		}
	}
	
	if ($list) {
		my $lid = $list =~ /^\d+$/ ? $list : $list->id;
		foreach my $l ( $self->lists() ) {
			if ( $lid == $l->id ) {
				return 1;
			}
		}		
	}
	
	if ($experiment) {
		my $eid = $experiment =~ /^\d+$/ ? $experiment : $experiment->id;
		foreach my $e ( $self->experiments(include_deleted => 1) ) {
			if ( $eid == $e->id ) {
				return 1;
			}
		}		
	}	
	
	return 0;
}

################################################ subroutine header begin ##

=head2 is_owner

 Usage     : $self->(dsg=>$dsg)
 Purpose   : checks to see if a user is the owner of a genome or dataset or list
 Returns   : 1/0
 Argument  : dsg=>genome object
             ds=>dataset object
             list=>list object
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub is_owner {
	my $self = shift;
	my %opts = @_;
	$opts{role} = 'owner';
	return is_role($self, %opts);
}

################################################ subroutine header begin ##

=head2 is_editor

 Usage     : $self->(dsg=>$dsg)
 Purpose   : checks to see if a user is the editor of a genome or dataset or list
 Returns   : 1/0
 Argument  : dsg=>genome object
             ds=>dataset object
             list=>list object
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub is_editor {
	my $self = shift;
	my %opts = @_;
	$opts{role} = 'editor';
	return is_role($self, %opts);
}

################################################ subroutine header begin ##

=head2 is_owner_editor

 Usage     : $self->(dsg=>$dsg)
 Purpose   : checks to see if a user is the owner/editor of a genome or dataset or list
 Returns   : 1/0
 Argument  : dsg=>genome object
             ds=>dataset object
             list=>list object
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub is_owner_editor {
	my $self = shift;
	my %opts = @_;

	$opts{role} = 'editor';
	return 1 if (is_role($self, %opts));

	$opts{role} = 'owner';
	return is_role($self, %opts);
}

################################################ subroutine header begin ##

=head2 is_reader

 Usage     : $self->(dsg=>$dsg)
 Purpose   : checks to see if a user is the reader of a genome or dataset or list
 Returns   : 1/0
 Argument  : dsg=>genome object
             ds=>dataset object
             list=>list object
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub is_reader {
	my $self = shift;
	my %opts = @_;
	$opts{role} = 'reader';
	return is_role($self, %opts);
}

################################################ subroutine header begin ##

=head2 is_role

 Usage     : $self->(dsg=>$dsg)
 Purpose   : checks to see if a user has given role in relation to given group/genome/dataset/list/experiment
 Returns   : 1/0
 Argument  : role => "owner","editor","reader" 
 			 group => user group object
             dsg => genome object
             ds => dataset object
             list => list object
             experiment => experiment object
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub is_role {
	my $self = shift;
	my %opts = @_;
	my $role = $opts{role};
	my $group = $opts{group};
	my $dsg  = $opts{dsg}; #FIXME rename to 'genome'
	my $ds   = $opts{ds};
	my $list = $opts{list};
	my $experiment = $opts{experiment};
	return 0 unless $self->id;
	return 0 unless  $dsg || $ds || $list || $experiment || $group;

	if ($group) {
		my $ugid = $group =~ /^\d+$/ ? $group : $group->id;
		foreach my $conn ($self->user_connectors) {
			return unless $conn->is_child_group;
			return 1 if ($conn->child_id == $ugid and $conn->role->name =~ /$role/i);
		}
	}

	if ($list) {
		my $lid = $list =~ /^\d+$/ ? $list : $list->id;
		foreach my $conn ($self->user_connectors) {
			next unless $conn->is_child_list;
			return 1 if ($conn->child_id == $lid and $conn->role->name =~ /$role/i);
		}
	}

	if ($dsg) { # genome
		my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
		foreach my $conn ($self->user_connectors) {
			next unless $conn->is_child_genome;
			return 1 if ($conn->child_id == $dsgid and $conn->role->name =~ /$role/i);
		}
	}

	if ($ds) { # dataset
		my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
		foreach my $conn ($self->user_connectors) {
			next unless $conn->is_child_genome;
			foreach ($conn->child->datasets) {
				return 1 if ($_->id == $dsid and $conn->role->name =~ /$role/i);
			}
		}
	}

	if ($experiment) {
		my $eid = $experiment =~ /^\d+$/ ? $experiment : $experiment->id;
		foreach my $conn ($self->user_connectors) {
			next unless $conn->is_child_experiment;
			return 1 if ($conn->child_id == $eid and $conn->role->name =~ /$role/i);
		}
	}

	# if ($group) {
	# 	return $group->role->name =~ /$role/i;
	# }

	# if ($dsg) {
	# 	my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
	# 	foreach my $group ( $self->groups ) {
	# 		next unless $group->role->name =~ /$role/i or $group->creator_user_id == $self->id;
	# 		foreach my $genome ( $group->genomes ) {
	# 			return 1 if $genome->id == $dsgid;
	# 		}
	# 	}
	# }

	# if ($ds) {
	# 	my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
	# 	foreach my $group ( $self->groups ) {
	# 		next unless $group->role->name =~ /$role/i or $group->creator_user_id == $self->id;
	# 		foreach my $ds ( $group->datasets ) {
	# 			return 1 if $ds->id == $dsid;
	# 		}
	# 	}
	# }
	
	# if ($list) {
	# 	my $lid = $list =~ /^\d+$/ ? $list : $list->id;
	# 	foreach my $group ( $self->groups ) {
	# 		next unless $group->role->name =~ /$role/i or $group->creator_user_id == $self->id;
	# 		foreach my $l ( $group->lists ) {
	# 			return 1 if $l->id == $lid;
	# 		}
	# 	}
	# }	
	
	# if ($experiment) {
	# 	my $eid = $experiment =~ /^\d+$/ ? $experiment : $experiment->id;
	# 	foreach my $group ( $self->groups ) {
	# 		next unless $group->role->name =~ /$role/i or $group->creator_user_id == $self->id;
	# 		foreach my $e ( $group->experiments ) {
	# 			return 1 if $e->id == $eid;
	# 		}
	# 	}
	# }	

	return 0;
}

################################################ subroutine header begin ##

=head2 datasets

 Usage     : 
 Purpose   : Returns the set of datasets a user has access to
 Returns   : Array of datasets
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub datasets {
	my $self = shift;
	return unless $self->id;
	
	my %datasets;
	foreach my $ug ( $self->groups ) {
		map { $datasets{ $_->id } = $_ } $ug->datasets;
	}
	return wantarray ? values %datasets : [ values %datasets ];		
}

################################################ subroutine header begin ##

=head2 restricted_datasets

 Usage     : 
 Purpose   : Returns the set of restricted datasets a user has access to
 Returns   : Array of datasets
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub restricted_datasets {
	my $self = shift;
	return unless $self->id;
	
	my %datasets;
	foreach my $ug ( $self->groups ) {
		map { $datasets{ $_->id } = $_ } $ug->restricted_datasets;
	}
	return wantarray ? values %datasets : [ values %datasets ];		
}

################################################ subroutine header begin ##

=head2 lists

 Usage     : $self->lists
 Purpose   : shows the lists to which user has access
 Returns   : wantarray of list objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub lists {
	my $self = shift;
	my %opts = @_;
	
	my %lists;
	foreach my $ug ( $self->groups ) { 
		#map { $lists{ $_->id } = $_ } $ug->lists; # FIXME will go away with new user_connector
		foreach my $uc ( $ug->child_connectors ) {
			if ($uc->is_child_list) {
				$lists{ $uc->child_id } = $uc->child;
			}
		}
	}
	foreach my $uc ( $self->user_connectors ) {
		if ($uc->is_child_list) {
			$lists{ $uc->child_id } = $uc->child;
		}
	}

	return wantarray ? values %lists : [ values %lists ];
}

################################################ subroutine header begin ##

=head2 experiments

 Usage     : $self->experiments
 Purpose   : Return set of experiments to which user has access
 Returns   : wantarray of experiment objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub experiments {
	my $self = shift;
	my %opts = @_;
	my $include_deleted = $opts{include_deleted};
	
	my %experiments;
	foreach my $ug ( $self->groups ) {
		# map { # FIXME will go away with new user_connector
		# 	$experiments{ $_->id } = $_ if (!$_->deleted || $include_deleted)
		# } $ug->experiments;
		
		foreach my $uc ( $ug->child_connectors ) {
			if ($uc->is_child_experiment) {
				my $experiment = $uc->child;
				next if ($experiment->deleted && not $include_deleted);
				$experiments{ $uc->child_id } = $experiment;
			}
			elsif ($uc->is_child_list) {
				my $list = $uc->child;
				map { $experiments{ $_->id } = $_ } $list->experiments( include_deleted => $include_deleted );
			}
		}
	}
	foreach my $uc ( $self->user_connectors ) {
		if ($uc->is_child_experiment) {
			my $experiment = $uc->child;
			next if ($experiment->deleted && not $include_deleted);
			$experiments{ $uc->child_id } = $experiment;
		}
		elsif ($uc->is_child_list) {
			my $list = $uc->child;
			map { $experiments{ $_->id } = $_ } $list->experiments( include_deleted => $include_deleted );
		}
	}

	return wantarray ? values %experiments : [ values %experiments ];	
}

################################################ subroutine header begin ##

=head2 restricted_experiments

 Usage     : $self->restricted_experiments
 Purpose   : Return set of restricted experiments to which user has access
 Returns   : wantarray of experiment objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub restricted_experiments {
	my $self = shift;
	my %opts = @_;
	my $include_deleted = $opts{include_deleted};
	
	my %experiments;
	foreach my $ug ( $self->groups ) {
		map { 
			$experiments{ $_->id } = $_ if (!$_->deleted || $include_deleted)
		} $ug->restricted_experiments;
	}
	return wantarray ? values %experiments : [ values %experiments ];	
}

################################################ subroutine header begin ##

=head2 genomes

 Usage     : $self->genomes
 Purpose   : Get list of genomes to which user has access
 Returns   : wantarray of genome objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub genomes {
	my $self = shift;
	my %opts = @_;
	my $include_deleted = $opts{include_deleted};
	my $ids = ($opts{ids}); #return genome ids only

	my %genomes;
	foreach my $ug ( $self->groups ) {
		# map { # FIXME will go away with new user_connector
		# 	$genomes{ $_->id } = $_ if (!$_->deleted || $include_deleted)
		# } $ug->genomes;
		
		foreach my $uc ( $ug->child_connectors ) {
			if ($uc->is_child_genome) {
###
			  if ($ids)
			    {
			      $genomes{$uc->child_id}=1;
			    }
			  else
			    {
			      my $genome = $uc->child;
			      next if ($genome->deleted && not $include_deleted);
			      $genomes{ $uc->child_id } = $genome;
			    }
##
			}
			elsif ($uc->is_child_list) {
				my $list = $uc->child;
				if ($ids)
				  {
				    map { $genomes{ $_ } = 1 } $list->genomes( include_deleted => $include_deleted, ids=>$ids );
				  }
				else
				  {
				    map { $genomes{ $_->id } = $_ } $list->genomes;
				  }
			}
		}
	}
	foreach my $uc ( $self->user_connectors ) {
		if ($uc->is_child_genome) {
###
			  if ($ids)
			    {
			      $genomes{$uc->child_id}=1;
			    }
			  else
			    {
			      my $genome = $uc->child;
			      next if ($genome->deleted && not $include_deleted);
			      $genomes{ $uc->child_id } = $genome;
			    }
##
		}
		elsif ($uc->is_child_list) {
			my $list = $uc->child;
			if ($ids)
			  {
			    map { $genomes{ $_ } = 1 } $list->genomes( include_deleted => $include_deleted, ids=>$ids );
			  }
			else
			  {
			    map { $genomes{ $_->id } = $_ } $list->genomes( include_deleted => $include_deleted, ids=>$ids );
			  }
		}
	}
	my @items;
	if ($ids) {@items = keys %genomes;}
	else {@items = values %genomes;}
	return wantarray ? @items : [ @items ];
}

################################################ subroutine header begin ##

=head2 restricted_genomes

 Usage     : 
 Purpose   : Returns the set of restricted genomes a user has access to
 Returns   : Array of genomes
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub restricted_genomes {
	my $self = shift;
	return unless $self->id;
	my %opts = @_;
	my $include_deleted = $opts{include_deleted};
	
	my %genomes;
	foreach my $ug ( $self->groups ) {
		map { 
			$genomes{ $_->id } = $_ if (!$_->deleted || $include_deleted)
		} $ug->restricted_genomes;
	}
	return wantarray ? values %genomes : [ values %genomes ];	
}

################################################ subroutine header begin ##

=head2 features

 Usage     : $self->features
 Purpose   : shows the features to which user has access
 Returns   : wantarray of features objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub features {
	my $self = shift;
#	my %opts = @_;
	
	my %features;
	foreach my $ug ( $self->groups ) {
		map { $features{ $_->id } = $_ } $ug->features;
	}
	return wantarray ? values %features : [ values %features ];
}

################################################ subroutine header begin ##

=head2 history

 Usage     : $self->history
 Purpose   : get the user's history
 Returns   : wantarray or count of history objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub history {
	my $self = shift;
	my %opts = @_;
	my $count = $opts{count}; #return count;

	if ($count) {
	    return $self->logs->count();
	}

	my @history = $self->logs;
	return wantarray ? @history : \@history;
}

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : generate a string of information about the user 
 Returns   : a string
 Argument  : None
 Throws    : None
 Comments  : uses name, description and role

=cut

################################################## subroutine header end ##

sub info {
	my $self = shift;
	return $self->user_name if (not $self->first_name and not $self->last_name); 
	return $self->first_name . ' ' . $self->last_name . ' (' . $self->user_name . ')';
}

1;

=head1 BUGS

=head1 SUPPORT

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

=cut

