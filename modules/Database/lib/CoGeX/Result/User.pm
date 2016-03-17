package CoGeX::Result::User;

use strict;
use warnings;

use base 'DBIx::Class::Core';
use base 'Class::Accessor';
use Data::Dumper;

=head1 NAME

CoGeX::User

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user> table in the CoGe database.

=head1 USAGE

  use CoGeX;

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen
 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("coge.user");
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
__PACKAGE__->belongs_to( 'image' => 'CoGeX::Result::Image', 'image_id');
__PACKAGE__->has_many( 'sessions' => "CoGeX::Result::UserSession", 'user_id' );
__PACKAGE__->has_many( 'works' => "CoGeX::Result::Work", 'user_id' );
__PACKAGE__->has_many( 'workflows' => "CoGeX::Result::Workflow", 'user_id' );
__PACKAGE__->has_many( 'logs' => "CoGeX::Result::Log", 'user_id' );
__PACKAGE__->has_many( 'jobs' => "CoGeX::Result::Job", 'user_id' );
__PACKAGE__->has_many( # all children (groups/genomes/experiments/lists)
	'child_connectors' => "CoGeX::Result::UserConnector",
	{ "foreign.parent_id" => "self.user_id" },
	{ where => { parent_type => $node_types->{user} } } );
__PACKAGE__->has_many( # child groups
	'group_connectors' => "CoGeX::Result::UserConnector",
	{ "foreign.parent_id" => "self.user_id" },
	{ where => [ -and => [ parent_type => $node_types->{user}, child_type => $node_types->{group} ] ] } );
__PACKAGE__->has_many( # child genomes
	'genome_connectors' => "CoGeX::Result::UserConnector",
	{ "foreign.parent_id" => "self.user_id" },
	{ where => [ -and => [ parent_type => $node_types->{user}, child_type => $node_types->{genome} ] ] } );
__PACKAGE__->has_many( # child experiments
	'experiment_connectors' => "CoGeX::Result::UserConnector",
	{ "foreign.parent_id" => "self.user_id" },
	{ where => [ -and => [ parent_type => $node_types->{user}, child_type => $node_types->{experiment} ] ] } );
__PACKAGE__->has_many( # child lists
	'list_connectors' => "CoGeX::Result::UserConnector",
	{ "foreign.parent_id" => "self.user_id" },
	{ where => [ -and => [ parent_type => $node_types->{user}, child_type => $node_types->{list} ] ] } );

__PACKAGE__->mk_accessors(qw(_genome_ids _experiment_ids));
#_genome_ids is a hash_ref of the genome_ids that a user has access to

sub item_type {
    return $node_types->{user};   
}

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
	return unless $self->id; # ignore public user
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

	my @groups;
	foreach my $conn ($self->group_connectors) {
		push @groups, $conn->child_group;
	}

	#print STDERR join(',', map {$_->id } @groups) . "\n";
	return wantarray ? @groups : \@groups;
}

sub groups {
	my $self = shift;
	return unless $self->id; # ignore public user
	return $self->user_groups(@_);
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
	return unless $self->id; # ignore public user

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
	return unless $self->id; # ignore public user
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
	return shift->admin;
}

sub is_public {
	my $self = shift;
	return !$self->id;
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
	my ($self, $list) = @_;
	return 0 unless $list;
	return 1 if (not $list->restricted or $self->is_admin); # public list or superuser
	return 0 if ($self->is_public); # deny public user for restricted list

	my $lid = $list->id;
	foreach ($self->lists) {
		return 1 if ($_->id == $lid);
	}

	return 0;
}

#new version caches results
sub has_access_to_genome {
	my ($self, $genome) = @_;
	return 0 unless $genome;
	return 1 if (not $genome->restricted or $self->is_admin); # public genome or superuser
	return 0 if ($self->is_public); # deny public user for restricted genome

	my $gid = $genome->id;
	unless ( $self->_genome_ids)  {
	    $self->_genome_ids({ map{$_->id=>1} $self->genomes(include_deleted => 1) });
	}
	my $ids = $self->_genome_ids();
	return $ids->{$gid};
}

#old version
#sub has_access_to_genome {
#	my $self = shift;
#	return unless $self->id; # ignore public user
#	my $genome  = shift;
#	my $gid = $genome =~ /^\d+$/ ? $genome : $genome->id;
#	return $self->is_admin || $self->child_connector(id => $gid, type => 'genome');
#}

# new version caches results
sub has_access_to_experiment {
	my ($self, $experiment) = @_;
	return 0 unless $experiment;
	return 1 if (not $experiment->restricted or $self->is_admin); # public experiment or superuser
	return 0 if ($self->is_public); # deny public user for restricted experiment

	unless ( $self->_experiment_ids) {
	    $self->_experiment_ids({ map{$_->id=>1} $self->experiments(include_deleted => 1) });
	}
	my $ids = $self->_experiment_ids();
	return $ids->{$experiment->id};
}

#old version
#sub has_access_to_experiment {
#	my $self = shift;
#	return unless $self->id; # ignore public user
#	my $experiment  = shift;
#	my $eid = $experiment =~ /^\d+$/ ? $experiment : $experiment->id;
#	return $self->is_admin || $self->child_connector(id => $eid, type => 'experiment');
#}

sub has_access_to_dataset {
	my ($self, $ds) = @_;
	return 0 unless $ds;
	return 1 if ($self->is_admin);#if (not $ds->restricted or $self->is_admin); # public dataset or superuser
	return 0 if ($self->is_public); # deny public user for restricted dataset

	my $dsid = $ds->id;
	foreach my $genome ($self->genomes(include_deleted => 1)) {
		foreach my $dataset ($genome->datasets) {
			return 1 if ($dataset->id == $dsid);
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
	$opts{role_id} = 2;#'owner';
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
	return unless $self->id; # ignore public user
	my %opts = @_;
	$opts{role_id} = 3;#'editor';
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
	return unless $self->id; # ignore public user
	my %opts = @_;

	$opts{role_id} = 3; #'editor';
	return 1 if (is_role($self, %opts));

	$opts{role_id} = 2; #'owner';
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
	return unless $self->id; # ignore public user
	my %opts = @_;
	$opts{role_id} = 4;#'reader';
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
	return 0 unless $self->id; # ignore public user
	my %opts = @_;
	#my $role = $opts{role};       # mdb removed 2/28/14
	my $role_id = $opts{role_id};  # mdb added 2/28/24
	my $group = $opts{group};
	my $dsg  = $opts{dsg}; #FIXME rename to 'genome'
	my $ds   = $opts{ds};
	my $list = $opts{list};
	my $experiment = $opts{experiment};
	return 0 unless $self->id; # not logged in
	return 0 unless $dsg || $ds || $list || $experiment || $group;

    if ($dsg) { # genome
        my $dsgid = $dsg =~ /^\d+$/ ? $dsg : $dsg->id;
        my $conn = $self->child_connector(id=>$dsgid, type=>'genome');
        #return 1 if ($conn && $conn->role->name =~ /$role/i); # mdb removed 2/28/14
        return 1 if ($conn && $conn->role_id == $role_id); # mdb added 2/28/14 for performance
    }

	if ($group) {
		my $ugid = $group =~ /^\d+$/ ? $group : $group->id;
		foreach my $conn ($self->group_connectors({child_id => $ugid})) {
			#return 1 if ($conn->role->name =~ /$role/i); # mdb removed 2/28/14
			return 1 if ($conn->role_id == $role_id); # mdb added 2/28/14 for performance
		}
	}

	if ($list) {
		my $lid = $list =~ /^\d+$/ ? $list : $list->id;
		my $conn = $self->child_connector(id=>$lid, type=>'list');
		#return 1 if ($conn && $conn->role->name =~ /$role/i); # mdb removed 2/28/14
		return 1 if ($conn && $conn->role_id == $role_id); # mdb added 2/28/14 for performance
	}

	if ($ds) { # dataset
		my $dsid = $ds =~ /^\d+$/ ? $ds : $ds->id;
		foreach my $conn ($self->all_child_connectors(type=>'genome')) {
			foreach ($conn->child->datasets) {
				#return 1 if ($_->id == $dsid && $conn->role->name =~ /$role/i); # mdb removed 2/28/14
				return 1 if ($_->id == $dsid && $conn->role_id == $role_id); # mdb added 2/28/14 for performance
			}
		}
	}

	if ($experiment) {
		my $eid = $experiment =~ /^\d+$/ ? $experiment : $experiment->id;
		my $conn = $self->child_connector(id=>$eid, type=>'experiment');
		#return 1 if ($conn && $conn->role->name =~ /$role/i); # mdb removed 2/28/14
		return 1 if ($conn && $conn->role_id == $role_id); # mdb added 2/28/14 for performance
	}

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
	return unless $self->id; # ignore public user

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
	return unless $self->id; # ignore public user

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
	return unless $self->id; # ignore public user
	my %opts = @_;

	my %lists;
	foreach my $ug ( $self->groups ) {
		foreach my $uc ( $ug->list_connectors ) {
			$lists{ $uc->child_id } = $uc->child;
		}
	}
	foreach my $uc ( $self->list_connectors ) {
		$lists{ $uc->child_id } = $uc->child;
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
	return unless $self->id; # ignore public user
	my %opts = @_;
	my $include_deleted = $opts{include_deleted}; # optional flag to include deleted

	my @experiments;
	foreach ($self->children(type=>'experiment')) {
		push @experiments, $_ unless ($_->deleted and not $include_deleted);
	}
	return wantarray ? @experiments : \@experiments;
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
	return unless $self->id; # ignore public user
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
	return unless $self->id; # ignore public user
	my %opts = @_;
	my $include_deleted = $opts{include_deleted}; # optional flag to include deleted genomes

	my @genomes;
	foreach ($self->children(type=>'genome')) {
		push @genomes, $_ unless ($_->deleted and not $include_deleted);
	}
	return wantarray ? @genomes : \@genomes;
}

sub groups_with_access {
	my $self = shift;
	return unless $self->id; # ignore public user
	my $item = shift;

	my %groups;

	# Groups
	foreach my $conn ( $item->group_connectors ) {
		my $group = $conn->parent;
		$groups{$group->id} = $group;
	}

	# Lists
	foreach my $conn ( $item->list_connectors ) {
		my $list = $conn->parent_list;
		foreach ($list->group_connectors) {
			my $group = $_->parent;
			$groups{$group->id} = $group;
		}
	}

	return wantarray ? values %groups : [ values %groups ];
}

sub users_with_access {
	my $self = shift;
	return unless $self->id; # ignore public user
	my $item = shift;

	my %users;

	# Direct user connections
	foreach my $conn ( $item->user_connectors ) {
		$users{$conn->parent_id} = $conn->parent;
	}

	# Users in groups
	foreach my $group ( $self->groups_with_access($item) ) {
		map { $users{$_->id} = $_ } $group->users;
	}

	# Lists
	foreach my $conn ( $item->list_connectors ) {
		my $list = $conn->parent_list;
		foreach ($list->user_connectors) {
			my $user = $_->parent;
            next unless defined $user;

			$users{$user->id} = $user;
		}
		foreach ($list->group_connectors) {
			my $group = $_->parent;
			map { $users{$_->id} = $_ } $group->users;
		}
	}

	return wantarray ? values %users : [ values %users ];
}

# Only call for children of type genome/experiment, not group/list.
sub child_connector {
	my $self = shift;
	return unless $self->id; # ignore public user
	my %opts = @_;
	my $id = $opts{id};
	my $type = $opts{type};
	my $type_num = $node_types->{$type};

	# Scan user's items - assumes there is only one user connector at this level
	foreach ($self->child_connectors({child_id=>$id, child_type=>$type_num})) {
		return $_;
	}

	# Scan user's lists
	if ($type ne 'list') { # Don't traverse lists within a list
		foreach my $conn ($self->list_connectors) {
			my $list = $conn->child;
			foreach ($list->child_connectors({child_id=>$id, child_type=>$type_num})) {
				next if ($_->child_id != $id or $_->child_type != $type_num); # FIXME mdb tempfix 10/14/13 issue 232 -- why necessary with line above?
				next unless ($conn->role_id == 4); # mdb tempfix 10/16/13 issue 232 -- list can only grant read access
				return $conn;
			}
		}
	}

	# Scan user's groups
	foreach my $group ($self->groups) { #TODO move this coge into UserGroup.pm::genomes ...?
		# Scan group's items
		foreach ($group->child_connectors({child_id=>$id, child_type=>$type_num})) {
			return $_;
		}
		# Scan group's lists
		if ($type ne 'list') { # Don't traverse lists within a list
			foreach my $conn ($group->list_connectors) {
				my $list = $conn->child;
				foreach ($list->child_connectors({child_id=>$id, child_type=>$type_num})) {
					next if ($_->child_id != $id or $_->child_type != $type_num); # FIXME mdb tempfix 10/14/13 issue 232 -- why necessary with line above?
					next unless ($conn->role_id == 4); # mdb tempfix 10/16/13 issue 232 -- list can only grant read access
					return $conn;
				}
			}
		}
	}
}

sub all_child_connectors { #FIXME optimize by mimicking child_by_type_and_id, combine with child_connector
	my $self = shift;
	return unless $self->id; # ignore public user
	my %opts = @_;
	my $type = $opts{type};
	my $type_num = $node_types->{$type};

	my %connectors;

	# Scan user's items
	foreach ($self->child_connectors({child_type=>$type_num})) {
		$connectors{$_->id} = $_;
	}

	# Scan user's lists
	if ($type ne 'list') {
		foreach my $conn ($self->list_connectors) {
			my $list = $conn->child;
			foreach ($list->child_connectors({child_type=>$type_num})) {
				$connectors{$_->id} = $_;
			}
		}
	}

	# Scan user's groups
	foreach my $group ($self->groups) { #TODO move this coge into UserGroup.pm::genomes ...?
		# Scan group's items
		foreach ($group->child_connectors({child_type=>$type_num})) {
			$connectors{$_->id} = $_;
		}
		# Scan group's lists
		if ($type ne 'list') {
			foreach my $conn ($group->list_connectors) {
				my $list = $conn->child;
				foreach ($list->child_connectors({child_type=>$type_num})) {
					$connectors{$_->id} = $_;
				}
			}
		}
	}

	return wantarray ? values %connectors : [ values %connectors ];
}

sub children { #FIXME have this use child_by_type_and_id

	my $self = shift;
	return unless $self->id; # ignore public user
	my %opts = @_;
	my $type = $opts{type};
	my $type_num = $node_types->{$type};

	my %children;

	# Scan user's items
	foreach ($self->child_connectors({child_type=>$type_num})) {
		my $child = $_->child;
		$children{$child->id} = $child;
	}

	# Scan user's lists
	foreach my $conn ($self->list_connectors) {
		my $list = $conn->child;
		foreach ($list->child_connectors({child_type=>$type_num})) {
			my $child = $_->child;
			$children{$child->id} = $child;
		}
	}

	# Scan user's groups
	foreach ($self->group_connectors) { #TODO move this coge into UserGroup.pm::genomes ...?
		my $group = $_->child;

		# Scan group's items
		foreach ($group->child_connectors({child_type=>$type_num})) {
			my $child = $_->child;
			$children{$child->id} = $child;
		}
		# Scan group's lists
		foreach my $conn ($group->list_connectors) {
			my $list = $conn->child;
			foreach ($list->child_connectors({child_type=>$type_num})) {
				my $child = $_->child;
				$children{$child->id} = $child;
			}
		}
	}

	return wantarray ? values %children : [ values %children ];
}

sub children_by_type_and_id {
	my $self = shift;
	return unless $self->id; # ignore public user
	my %opts = @_;

	# use Time::HiRes qw ( time );
	# my $start_time = time;

	my %children;

	foreach my $c ($self->child_connectors) {
		my $child = $c->child;
		$children{$c->child_type}{$c->child_id} = $child;

		if ($c->child_type == $node_types->{list}) {
			foreach my $c ($child->child_connectors) {
				my $child = $c->child;
				$children{$c->child_type}{$c->child_id} = $child;
			}
		}
		elsif ($c->child_type == $node_types->{group}) {
			foreach my $c ($child->child_connectors) {
				my $child = $c->child;
				$children{$c->child_type}{$c->child_id} = $child;

				if ($c->child_type == $node_types->{list}) {
					foreach my $c ($child->child_connectors) {
						my $child = $c->child;
						$children{$c->child_type}{$c->child_id} = $child;
					}
				}
			}
		}
	}

#	print STDERR "children_by_type_and_id: time=" . ((time - $start_time)*1000) . "\n";

	return \%children;
}

sub children_by_type_role_id {
	my $self = shift;
	return unless $self->id; # ignore public user

	#use Time::HiRes qw ( time ); my $start_time = time;
	my (%children, %roles);

	foreach my $c ($self->child_connectors) {
		my $child = $c->child;
		$children{$c->child_type}{$c->child_id} = $child;
		$roles{$c->role_id}{$c->child_id} = 1;

		if ($c->child_type == $node_types->{list}) {
			foreach my $lc ($child->child_connectors) {
				my $child = $lc->child;
				$children{$lc->child_type}{$lc->child_id} = $child;
				$roles{$c->role_id}{$lc->child_id} = 1;
			}
		}
		elsif ($c->child_type == $node_types->{group}) {
			foreach my $c ($child->child_connectors) {
				my $child = $c->child;
				$children{$c->child_type}{$c->child_id} = $child;
				$roles{$c->role_id}{$c->child_id} = 1;

				if ($c->child_type == $node_types->{list}) {
					foreach my $lc ($child->child_connectors) {
						my $child = $lc->child;
						$children{$lc->child_type}{$lc->child_id} = $child;
						$roles{$c->role_id}{$lc->child_id} = 1;
					}
				}
			}
		}
	}
	#print STDERR "children_by_type_and_id: time=" . ((time - $start_time)*1000) . "\n";

	return (\%children, \%roles);
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
	return unless $self->id; # ignore public user
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
	return unless $self->id; # ignore public user
#	my %opts = @_;

	my %features;
	foreach my $group ( $self->groups ) {
		map { $features{ $_->id } = $_ } $group->features;
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
	return unless $self->id; # ignore public user
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
	return unless $self->id; # ignore public user
	my $name = $self->user_name . ', id' . $self->id;
	return $name if (not $self->first_name and not $self->last_name);
	return $self->first_name . ' ' . $self->last_name . ' (' . $name . ')';
}

1;
