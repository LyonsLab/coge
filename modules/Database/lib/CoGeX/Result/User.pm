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

__PACKAGE__->table("user");
__PACKAGE__->add_columns(
	"user_id",
	{
		data_type     => "INT",
		default_value => undef,
		is_nullable   => 0,
		size          => 10
	},
	"user_name",
	{
		data_type     => "VARCHAR",
		default_value => "",
		is_nullable   => 0,
		size          => 10
	},
	"first_name",
	{
		data_type     => "VARCHAR",
		default_value => "",
		is_nullable   => 0,
		size          => 10
	},
	"last_name",
	{
		data_type     => "VARCHAR",
		default_value => "",
		is_nullable   => 0,
		size          => 10
	},
	"email",
	{
		data_type     => "VARCHAR",
		default_value => undef,
		is_nullable   => 1,
		size          => 50,
	},
	"description",
	{
		data_type     => "VARCHAR",
		default_value => undef,
		is_nullable   => 1,
		size          => 255,
	},
);
__PACKAGE__->set_primary_key("user_id");
__PACKAGE__->has_many( 'sessions'  => "CoGeX::Result::UserSession", 'user_id' );
__PACKAGE__->has_many( 'works'     => "CoGeX::Result::Work",        'user_id' );
__PACKAGE__->has_many( 'workflows' => "CoGeX::Result::Workflow",    'user_id' );
__PACKAGE__->has_many(
	'user_group_connectors' => "CoGeX::Result::UserGroupConnector",
	'user_id'
);

#__PACKAGE__->has_many('lists'=>"CoGeX::Result::List",'user_id');

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
	return unless $self->id;

	my @user_groups = ();

	foreach my $user_group_connector ( $self->user_group_connectors() ) {
		push( @user_groups, $user_group_connector->user_group() )
		  if $user_group_connector->user_group();
	}

	return wantarray ? @user_groups : \@user_groups;
}

################################################ subroutine header begin ##

=head2 groups

 Usage     : alias for $self->user_groups
 Purpose   : 
 Returns   : array or arrayref of user_group objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub groups {
	return shift->user_groups(@_);
}

################################################ subroutine header begin ##

=head2 owner_group

 Usage     : return user's locked owner group
 Purpose   : 
 Returns   : user_group object
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub owner_group {
	my $self = shift;
	foreach my $group ( $self->groups ) {
		return $group if ($group->locked && $group->role->name =~ /owner/i);
	}
	return;
}

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
	my $self = shift;
	foreach my $group ( $self->groups ) {
		return 1 if $group->role->name =~ /admin/i;
	}
	return 0;
}

################################################ subroutine header begin ##

=head2 add_to_group

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : Group Name
 Throws    : None
 Comments  : 



=cut

################################################## subroutine header end ##

sub add_to_group() {
	my $self = shift;

	my %opts       = @_;
	my $coge       = $opts{coge};
	my $group_name = $opts{group_name};

	my $group = $coge->ResultSet('UserGroup')->find( { name => $group_name } );

	if ($group) {

		# FIXME
	}
}

################################################ subroutine header begin ##

=head2 private_genomes

 Usage     : 
 Purpose   : Returns the set of private genomes a user has access to
 Returns   : Array of dataset_groups
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub private_genomes {

	my $self = shift;
	return unless $self->id;
	my @private_genomes = ();
	foreach my $group ( $self->user_groups() ) {
		push( @private_genomes, $group->private_genomes() );
	}

	return wantarray ? @private_genomes : \@private_genomes;
}

################################################ subroutine header begin ##

=head2 has_access_to_genome

 Usage     : 
 Purpose   : checks to see if a user has access to a genome (dataset_group)
 Returns   : 1/0
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub has_access_to_genome { # FIXME mdb 8/10/12 possibly obsolete
	my $self = shift @_;
	my $dsg  = shift;
	return 0 unless $dsg;
	foreach my $group ( $self->groups ) {
		return 1 if $group->role->name eq "Admin";
	}
	my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
	foreach my $genome ( $self->private_genomes() ) {
		if ( $dsgid == $genome->id ) {
			return 1;
		}
	}
	return 0;
}

################################################ subroutine header begin ##

=head2 private_datasets

 Usage     : 
 Purpose   : Returns the set of private datasets a user has access to
 Returns   : Array of datasets
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub private_datasets {

	my $self = shift;
	return unless $self->id;
	my @private_datasets = ();
	foreach my $group ( $self->user_groups() ) {
		push( @private_datasets, $group->private_datasets() );
	}

	return wantarray ? @private_datasets : \@private_datasets;
}

################################################ subroutine header begin ##

=head2 has_access_to_dataset

 Usage     : 
 Purpose   : checks to see if a user has access to a dataset
 Returns   : 1/0
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub has_access_to_dataset { # FIXME mdb 8/10/12 possibly obsolete
	my $self = shift @_;
	my $ds   = shift;
	return 0 unless $ds;
	foreach my $group ( $self->groups ) {
		return 1 if $group->role->name eq "Admin";
	}
	my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
	foreach my $ds ( $self->private_datasets() ) {
		if ( $dsid == $ds->id ) {
			return 1;
		}
	}
	return 0;
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

sub has_access { # FIXME mdb 8/10/12 possibly obsolete
	my $self = shift;
	my %opts = @_;
	my $dsg  = $opts{dsg};
	my $ds   = $opts{ds};
	my $list = $opts{list};
	my $experiment = $opts{experiment};
	return 0 unless $dsg || $ds || $list || $experiment;
	
	foreach my $group ( $self->groups ) {
		return 1 if $group->role->name eq "Admin";
	}
	
	if ($ds) {
		my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
		foreach my $ds ( $self->private_datasets() ) {
			if ( $dsid == $ds->id ) {
				return 1;
			}
		}
	}
	
	if ($dsg) {
		my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
		foreach my $genome ( $self->private_genomes() ) {
			if ( $dsgid == $genome->id ) {
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
		foreach my $e ( $self->experiments() ) {
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
	my $is_editor = is_role($self, %opts);
	return $is_editor if ($is_editor);

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
 Purpose   : checks to see if a user has given role in relation to given genome/dataset/list/experiment
 Returns   : 1/0
 Argument  : role => "owner","editor","reader" 
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
	my $dsg  = $opts{dsg}; #FIXME rename to 'genome'
	my $ds   = $opts{ds};
	my $list = $opts{list};
	my $experiment = $opts{experiment};
	return 0 unless $dsg || $ds || $list || $experiment;

	if ($dsg) {
		my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
		foreach my $group ( $self->groups ) {
			next unless $group->role->name =~ /$role/i;
			foreach my $genome ( $group->genomes ) {
				return 1 if $genome->id == $dsgid;
			}
		}
	}

	if ($ds) {
		my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
		foreach my $group ( $self->groups ) {
			next unless $group->role->name =~ /$role/i;
			foreach my $ds ( $group->datasets ) {
				return 1 if $ds->id == $dsid;
			}
		}
	}
	
	if ($list) {
		my $lid = $list =~ /^\d+$/ ? $list : $list->id;
		foreach my $group ( $self->groups ) {
			next unless $group->role->name =~ /$role/i;
			foreach my $l ( $group->lists ) {
				return 1 if $l->id == $lid;
			}
		}
	}	
	
	if ($experiment) {
		my $eid = $experiment =~ /^\d+$/ ? $experiment : $experiment->id;
		foreach my $group ( $self->groups ) {
			next unless $group->role->name =~ /$role/i;
			foreach my $e ( $group->experiments ) {
				return 1 if $e->id == $eid;
			}
		}
	}	

	return 0;
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
		map { $lists{ $_->id } = $_ } $ug->lists;
	}
	return wantarray ? values %lists : [ values %lists ];
}

################################################ subroutine header begin ##

=head2 experiments

 Usage     : $self->experiments
 Purpose   : shows the experiments to which user has access
 Returns   : wantarray of experiment objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub experiments {
	my $self = shift;
	my %opts = @_;
	my $restricted = $opts{restricted};
	
	my %experiments;
	foreach my $ug ( $self->groups ) {
		foreach my $list ( $ug->lists) {
			map { $experiments{ $_->id } = $_ } $list->experiments(restricted => $restricted);
		}
	}
	return wantarray ? values %experiments : [ values %experiments ];
}

################################################ subroutine header begin ##

=head2 genomes

 Usage     : $self->genomes
 Purpose   : shows the genomes to which user has access
 Returns   : wantarray of genome objects
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub genomes {
	my $self = shift;
	my %opts = @_;
	my %genomes;
	foreach my $ug ( $self->groups ) {
		foreach my $list ( $ug->lists) {
			map { $genomes{ $_->id } = $_ } $list->genomes;
		}
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
	my %opts = @_;
	my %features;
	foreach my $ug ( $self->groups ) {
		foreach my $list ( $ug->lists) {
			map { $features{ $_->id } = $_ } $list->features;
		}
	}
	return wantarray ? values %features : [ values %features ];
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
	my $info = $self->user_name;
	$info .= ": " . $self->first_name;
	$info .= " " . $self->last_name;
	return $info;
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

