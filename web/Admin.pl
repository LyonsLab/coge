#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON qw(encode_json);
use Data::Dumper;
use List::Compare;
use CoGeX;
use CoGeX::Result::User;
use Data::Dumper;
use URI::Escape::JavaScript qw(escape unescape);
no warnings 'redefine';

use vars
  qw($P $PAGE_NAME $USER $BASEFILE $coge $cogeweb %FUNCTION $FORM $MAX_SEARCH_RESULTS %ITEM_TYPE);

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init( cgi => $FORM );

$MAX_SEARCH_RESULTS = 400;

my $node_types = CoGeX::node_types();

#print STDERR $node_types->{user};

%FUNCTION = (
	search_organisms                => \&search_organisms,
	search_users                    => \&search_users,
	search_stuff                    => \&search_stuff,
	user_info                       => \&user_info,
	add_items_to_user_or_group      => \&add_items_to_user_or_group,
	remove_items_from_user_or_group => \&remove_items_from_user_or_group,
	get_share_dialog                => \&get_share_dialog,
	get_roles                       => \&get_roles,
	search_share                    => \&search_share,
	delete_genome                   => \&delete_genome,
	delete_list                     => \&delete_list,
	delete_experiment               => \&delete_experiment,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
	#print STDERR "HTML\n";
	my $html;
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );

	#$template->param( HELP => '/wiki/index.php?title=ADMIN' );
	$template->param( HELP => $P->{SERVER} );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( PAGE_TITLE => qq{Admin} );
	$template->param( TITLE      => "Admin" );
	$template->param( LOGO_PNG   => "CoGe.svg" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( BODY       => gen_body() );
	$template->param( ADJUST_BOX => 1 );
	$template->param( CAS_URL    => $P->{CAS_URL} || '' );
	$html .= $template->output;
}

sub gen_body {

	#print STDERR "BODY\n";
	# Hide this page if the user is not an Admin
	unless ( $USER->is_admin ) {
		my $template =
		  HTML::Template->new( filename => $P->{TMPLDIR} . "Admin.tmpl" );
		$template->param( ADMIN_ONLY => 1 );
		return $template->output;
	}

	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'Admin.tmpl' );
	$template->param( MAIN => 1 );

	#print STDERR $node_types->{user};
	#$template->param( ITEM_TYPE_USER => $node_types->{user} );
	return $template->output;
}

sub search_stuff {
	my %opts        = @_;
	my $search_term = $opts{search_term};
	my $timestamp   = $opts{timestamp};

	my @searchArray = split( ' ', $search_term );
	my @specialTerms;
	my @idList;
	my @results;

	#return unless $search_term;

	#Set up the necessary arrays of serch terms and special conditions
	for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
		if ( index( $searchArray[$i], ':' ) == -1 ) {
			$searchArray[$i] = { 'like', '%' . $searchArray[$i] . '%' };
		}
		else {
			my @splitTerm = split( ':', $searchArray[$i] );
			splice( @searchArray, $i, 1 );
			$i--;
			push @specialTerms,
			  { 'tag' => $splitTerm[0], 'term' => $splitTerm[1] };
		}
	}

	#say STDERR Dumper(\@searchArray);

	#Set the special conditions
	my $type = "none";    #Specifies a particular field to show
	my @restricted =
	  [ -or => [ restricted => 0, restricted => 1 ] ]; 
	  #Specifies either restricted(1) OR unrestriced(0) results. Default is all.
	my @deleted =
	  [ -or => [ deleted => 0, deleted => 1 ] ]; 
	  #Specifies either deleted(1) OR non-deleted(0) results. Default is all.
	for ( my $i = 0 ; $i < @specialTerms ; $i++ ) {
		if ( $specialTerms[$i]{tag} eq 'type' ) {
			$type = $specialTerms[$i]{term};
		}
		if ( $specialTerms[$i]{tag} eq 'restricted' ) {
			@restricted = [ restricted => $specialTerms[$i]{term} ];
			if ( $type eq "none" ) {
				$type = 'restricted'; 
				#Sets the "type" so that only fields relevant to restriction are shown. (i.e. there are no restricted users.)
			}
		}
		if (
			$specialTerms[$i]{tag} eq 'deleted'
			&& (   $specialTerms[$i]{term} eq '0'
				|| $specialTerms[$i]{term} eq '1' )
		  )
		{
			@deleted = [ deleted => $specialTerms[$i]{term} ];
			if ( $type eq "none" ) {
				$type = 'deleted'; 
				#Sets the "type" so that only fields relevant to deletion are shown. (i.e. there are no deleted users).
			}
		}

        #if ($specialTerms[$i]{tag} eq 'deleted' && $specialTerms[$i]{term} eq '*') {
        #        @deleted = [-or => [deleted => 0, deleted => 1]];
        #	if($type eq "none") {
        #                $type = 'deleted';
        #        }
        #}
	}

	# Perform organism search
	if (   $type eq 'none'
		|| $type eq 'organism'
		|| $type eq 'genome'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
		my @orgArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			push @orgArray,
			  [
				-or => [
					name        => $searchArray[$i],
					description => $searchArray[$i],
					organism_id => $searchArray[$i]
				]
			  ];
		}
		my @organisms =
		  $coge->resultset("Organism")->search( { -and => [ @orgArray, ], } );

		if ( $type eq 'none' || $type eq 'organism' ) {
			foreach ( sort { $a->name cmp $b->name } @organisms ) {
				push @results,
				  { 'type' => "organism", 'label' => $_->name, 'id' => $_->id };
			}
		}
		@idList = map { $_->id } @organisms;
	}

	# Perform user search
	if ( $type eq 'none' || $type eq 'user' ) {
		my @usrArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			push @usrArray,
			  [
				-or => [
					user_name  => $searchArray[$i],
					first_name => $searchArray[$i],
					user_id    => $searchArray[$i],
					first_name => $searchArray[$i]
				]
			  ];
		}
		my @users =
		  $coge->resultset("User")->search( { -and => [ @usrArray, ], } );

		foreach ( sort { $a->user_name cmp $b->user_name } @users ) {
			push @results,
			  { 'type' => "user", 'label' => $_->user_name, 'id' => $_->id };
		}
	}

	# Perform genome search (corresponding to Organism results)
	if (   $type eq 'none'
		|| $type eq 'genome'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
		my @genomes = $coge->resultset("Genome")->search(
			{
				-and => [
					organism_id => { -in => \@idList },
					@restricted,
					@deleted,
				],
			}
		);

		foreach ( sort { $a->id cmp $b->id } @genomes ) {
			push @results,
			  {
				'type'    => "genome",
				'label'   => $_->info,
				'id'      => $_->id,
				'deleted' => $_->deleted
			  };
		}

		# Perform direct genome search (by genome ID)
		my @genIDArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			push @genIDArray, [ -or => [ genome_id => $searchArray[$i] ] ];
		}
		my @genomeIDs =
		  $coge->resultset("Genome")
		  ->search( { -and => [ @genIDArray, @restricted, @deleted, ], } );

		foreach ( sort { $a->id cmp $b->id } @genomeIDs ) {
			push @results,
			  {
				'type'    => "genome",
				'label'   => $_->info,
				'id'      => $_->id,
				'deleted' => $_->deleted
			  };
		}
	}

	# Perform experiment search
	if (   $type eq 'none'
		|| $type eq 'experiment'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
		my @expArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			push @expArray,
			  [
				-or => [
					name          => $searchArray[$i],
					description   => $searchArray[$i],
					experiment_id => $searchArray[$i]
				]
			  ];
		}
		my @experiments =
		  $coge->resultset("Experiment")
		  ->search( { -and => [ @expArray, @restricted, @deleted, ], } );

		foreach ( sort { $a->name cmp $b->name } @experiments ) {
			push @results,
			  {
				'type'    => "experiment",
				'label'   => $_->name,
				'id'      => $_->id,
				'deleted' => $_->deleted
			  };
		}
	}

	# Perform notebook search
	if (   $type eq 'none'
		|| $type eq 'notebook'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
		my @noteArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			push @noteArray,
			  [
				-or => [
					name        => $searchArray[$i],
					description => $searchArray[$i],
					list_id     => $searchArray[$i]
				]
			  ];
		}
		my @notebooks =
		  $coge->resultset("List")
		  ->search( { -and => [ @noteArray, @restricted, @deleted, ], } );

		foreach ( sort { $a->name cmp $b->name } @notebooks ) {
			push @results,
			  {
				'type'    => "notebook",
				'label'   => $_->info,
				'id'      => $_->id,
				'deleted' => $_->deleted
			  };
		}
	}

	# Perform user group search
	if ( $type eq 'none' || $type eq 'usergroup' || $type eq 'deleted' ) {
		my @usrGArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			push @usrGArray,
			  [
				-or => [
					name          => $searchArray[$i],
					description   => $searchArray[$i],
					user_group_id => $searchArray[$i]
				]
			  ];
		}
		my @userGroup =
		  $coge->resultset("UserGroup")
		  ->search( { -and => [ @usrGArray, @deleted, ], } );

		foreach ( sort { $a->name cmp $b->name } @userGroup ) {
			push @results,
			  {
				'type'    => "user_group",
				'label'   => $_->name,
				'id'      => $_->id,
				'deleted' => $_->deleted
			  };
		}
	}

	#print STDERR "Successful search";
	return encode_json( { timestamp => $timestamp, items => \@results } );
}

sub user_info {

	#print STDERR "Hello? Is anyone there?";
	#print STDERR "@_\n";

	my %opts        = @_;
	my $search_term = $opts{search_term};
	my $search_type = $opts{search_type};
	my $timestamp   = $opts{timestamp};

	#print STDERR "$search_term\n";

	my $user;
	if ( $search_type eq "user" ) {
		$user = $coge->resultset("User")->find($search_term);
	}
	if ( $search_type eq "group" ) {
		$user = $coge->resultset("UserGroup")->find($search_term);
	}
	my @results;
	my @users;
	push( @users, $user );

	my $child;
	if ( $search_type eq "user" ) {
		foreach ( $user->child_connectors( { child_type => 6 } ) ) {
			$child = $_->child;
			push( @users, $child );

#push @results, { 'type' => "user_group", 'label' => $child->name, 'id' => $child->id, 'deleted' => $child->deleted};
		}
	}

	foreach my $currentUser (@users) {
		my @current;

		# Find notebooks
		foreach ( $currentUser->child_connectors( { child_type => 1 } ) ) {
			$child = $_->child;

#push @results, { 'type' => "notebook", 'label' => $child->name, 'id' => $child->id, 'info' => $child->info};
			push @current,
			  {
				'type'    => "notebook",
				'label'   => $child->info,
				'id'      => $child->id,
				'role'    => $_->role_id,
				'deleted' => $child->deleted
			  };
		}

		# Find genomes
		foreach ( $currentUser->child_connectors( { child_type => 2 } ) ) {
			$child = $_->child;
			push @current,
			  {
				'type'    => "genome",
				'label'   => $child->info,
				'id'      => $child->id,
				'role'    => $_->role_id,
				'deleted' => $child->deleted
			  };
		}

		# Find experiments
		foreach ( $currentUser->child_connectors( { child_type => 3 } ) ) {
			$child = $_->child;
			push @current,
			  {
				'type'    => "experiment",
				'label'   => $child->name,
				'id'      => $child->id,
				'info'    => $child->info,
				'role'    => $_->role_id,
				'deleted' => $child->deleted
			  };

#push @current, { 'type' => "experiment", 'label' => $child->info, 'id' => $child->id, 'role' => $_->role_id, 'deleted' => $child->deleted};
		}

		# Find users if searching a user group
		if ( $search_type eq "group" ) {
			foreach ( $user->users ) {
				push @current,
				  { 'type' => "user", 'label' => $_->name, 'id' => $_->id };
			}
		}

		push @results,
		  {
			'user'    => $currentUser->name,
			'user_id' => $currentUser->id,
			'result'  => \@current
		  };
	}

	#print STDERR Dumper(@results);
	return encode_json( { timestamp => $timestamp, items => \@results } );
}

sub add_items_to_user_or_group {
	my %opts        = @_;
	my $target_item = $opts{target_item};
	return unless $target_item;
	my $role_id = $opts{role_id};
	return unless $role_id;
	my $item_list = $opts{item_list};
	my @items = split( ',', $item_list );
	return unless @items;

	#print STDERR Dumper(\@items);
	#print STDERR "\n$target_item\n";
	#print STDERR "$role_id\n";

	my ( $target_id, $target_type ) = $target_item =~ /(\d+)\:(\d+)/;

	#print STDERR "$target_id $target_type\n";
	return unless ( $target_id and $target_type );

	# Verify that user has access to each item
	my @verified;
	foreach my $item (@items) {
		my ( $item_id, $item_type ) = $item =~ /content_(\d+)_(\d+)/;

		#print STDERR "$item_id $item_type\n";
		next unless ( $item_id and $item_type );

		# print STDERR "add_items_to_user_or_group $item_id $item_type\n";
		if ( $item_type == $node_types->{genome} ) {
			my $genome = $coge->resultset('Genome')->find($item_id);
			next
			  unless ( $USER->has_access_to_genome($genome)
				|| $USER->is_admin );
			push @verified, $item;
		}
		elsif ( $item_type == $node_types->{experiment} ) {
			my $experiment = $coge->resultset('Experiment')->find($item_id);
			next unless $USER->has_access_to_experiment($experiment);
			push @verified, $item;
		}
		elsif ( $item_type == $node_types->{notebook} ) {
			my $notebook = $coge->resultset('List')->find($item_id);
			next unless $USER->has_access_to_list($notebook);
			push @verified, $item;
		}
	}

	#print STDERR Dumper(\@verified);

	# Assign each item to user/group
	# print STDERR "add_items_to_user_or_group $target_id $target_type\n";
	#TODO verify that user can use specified role (for admin/owner roles)
	if ( $target_type == $node_types->{user} ) {
		my $user = $coge->resultset('User')->find($target_id);

		#print STDERR "$user\n";
		return unless $user;

		foreach (@verified) {
			my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

			#print STDERR "   user: $item_id $item_type\n";

			# Remove previous connection
			foreach (
				$coge->resultset('UserConnector')->search(
					{
						parent_id   => $target_id,
						parent_type => 5,            # FIXME hardcoded
						child_id    => $item_id,
						child_type  => $item_type
					}
				)
			  )
			{
				$_->delete;
			}

			# Add new connection
			my $conn = $coge->resultset('UserConnector')->create(
				{
					parent_id   => $target_id,
					parent_type => 5,            # FIXME hardcoded
					child_id    => $item_id,
					child_type  => $item_type,
					role_id     => $role_id
				}
			);

			#print STDERR "$conn\n";
			return unless $conn;
		}
	}
	elsif ( $target_type == $node_types->{group} ) {
		my $group = $coge->resultset('UserGroup')->find($target_id);

		return unless $group;

		foreach (@verified) {
			my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

			#print STDERR "   group: $item_id $item_type\n";
			#my $var = $group->role_id;
			#print STDERR "$target_id $var\n";

			my $conn = $coge->resultset('UserConnector')->find_or_create(
				{
					parent_id   => $target_id,
					parent_type => 6,                # FIXME hardcoded
					child_id    => $item_id,
					child_type  => $item_type,
					role_id     => $group->role_id
				}
			);
			return unless $conn;
		}
	}

	return get_share_dialog( item_list => $item_list );
}

sub remove_items_from_user_or_group {
	my %opts        = @_;
	my $target_item = $opts{target_item};

	#print STDERR "$target_item\n";
	return unless $target_item;
	my $item_list = $opts{item_list};
	my @items = split( ',', $item_list );

	#print STDERR Dumper(\@items);
	return unless @items;

	my ( $target_id, $target_type ) = $target_item =~ /(\d+)\:(\d+)/;

	#print STDERR "remove target $target_id $target_type\n";
	next unless ( $target_id and $target_type );

	foreach (@items) {
		my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

		#print STDERR "remove_item_from_user $item_id $item_type\n";
		next unless ( $item_id and $item_type );

		if ( $item_type == $node_types->{genome} ) {
			my $genome = $coge->resultset('Genome')->find($item_id);
			next unless ( $USER->has_access_to_genome($genome) );

			my $conn = $coge->resultset('UserConnector')->find(
				{
					parent_id   => $target_id,
					parent_type => $target_type,
					child_id    => $genome->id,
					child_type  => $node_types->{genome}
				}
			);
			return unless $conn;

			$conn->delete;
		}
		elsif ( $item_type == $node_types->{experiment} ) {
			my $experiment = $coge->resultset('Experiment')->find($item_id);
			next unless $USER->has_access_to_experiment($experiment);

			my $conn = $coge->resultset('UserConnector')->find(
				{
					parent_id   => $target_id,
					parent_type => $target_type,               #FIXME hardcoded
					child_id    => $experiment->id,
					child_type  => $node_types->{experiment}
				}
			);
			return unless $conn;

			$conn->delete;
		}
		elsif ( $item_type == $node_types->{notebook} ) {
			my $notebook = $coge->resultset('List')->find($item_id);
			next unless $USER->has_access_to_list($notebook);

			my $conn = $coge->resultset('UserConnector')->find(
				{
					parent_id   => $target_id,
					parent_type => $target_type,             #FIXME hardcoded
					child_id    => $notebook->id,
					child_type  => $node_types->{notebook}
				}
			);
			return unless $conn;

			$conn->delete;
		}
	}

	return get_share_dialog( item_list => $item_list );
}

sub get_share_dialog {    #FIXME this routine needs to be optimized
	                      #print STDERR "GetShare called.\n";
	my %opts      = @_;
	my $item_list = $opts{item_list};

	#print STDERR Dumper($item_list);
	my @items = split( ',', $item_list );
	return unless @items;

	my ( %userconn, %notebooks );
	my $isPublic   = 0;
	my $isEditable = 1;
	my $item_id;
	my $item_type;
	foreach (@items) {
		( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

		#print STDERR $item_id;
		#print STDERR "\n!";
		#print STDERR $item_type;
		#print STDERR "\n!!";
		next unless ( $item_id and $item_type );

		# print STDERR "get_share $item_id $item_type\n";
		if ( $item_type == $node_types->{genome} ) {
			my $genome = $coge->resultset('Genome')->find($item_id);
			next unless ( $USER->has_access_to_genome($genome) );
			map { $userconn{ $_->parent_id } = $_ }
			  ( $genome->user_connectors, $genome->group_connectors );
			map { $notebooks{ $_->id } = $_ } $genome->lists;
			$isPublic = 1 if ( not $genome->restricted );
			$isEditable = 0 if ( not $USER->is_owner_editor( dsg => $genome ) );
		}
		elsif ( $item_type == $node_types->{experiment} ) {
			my $experiment = $coge->resultset('Experiment')->find($item_id);
			next unless $USER->has_access_to_experiment($experiment);
			map { $userconn{ $_->id } = $_ }
			  ( $experiment->user_connectors, $experiment->group_connectors );
			map { $notebooks{ $_->id } = $_ } $experiment->lists;
			$isPublic = 1 if ( not $experiment->restricted );
			$isEditable = 0
			  if ( not $USER->is_owner_editor( experiment => $experiment ) );
		}
		elsif ( $item_type == $node_types->{notebook} ) {
			my $notebook = $coge->resultset('List')->find($item_id);
			next unless $USER->has_access_to_list($notebook);
			map { $userconn{ $_->id } = $_ }
			  ( $notebook->user_connectors, $notebook->group_connectors );
			$isPublic = 1 if ( not $notebook->restricted );
			$isEditable = 0
			  if ( not $USER->is_owner_editor( list => $notebook ) );
		}
	}

	$isEditable = 1 if ( $USER->is_admin );

	my ( %user_rows, %group_rows, %notebook_rows );
	foreach my $conn ( values %userconn ) {
		next
		  unless $conn->role
		; #EL added 10/21/14 so solve a problem for a user where the share dialog wouldn't appear for his genomes, and an fatal error was being thrown due to role->name with role being undefined.  Genomes with problems were IDs: 24576, 24721, 24518, 24515, 24564, 24566, 24568, 24562, 24571
		if ( $conn->is_parent_user ) {
			my $user = $conn->parent;

			#print STDERR $user;
			#print STDERR "\nYo\n";
			$user_rows{ $user->id } = {
				ITEM_ID        => $item_id,
				ITEM_TYPE      => $item_type,
				USER_ITEM      => $user->id . ':' . $conn->parent_type,
				USER_FULL_NAME => $user->display_name,
				USER_NAME      => $user->name,
				USER_ROLE      => $conn->role->name,
				USER_DELETE    => $isEditable
				  && (!$conn->role->is_owner
					|| $USER->is_admin )   # owner can't be removed unless admin
			};
		}
		elsif ( $conn->is_parent_group ) {
			my $group = $conn->parent;

			my @users = map {
				{
					GROUP_USER_FULL_NAME => $_->display_name,
					GROUP_USER_NAME      => $_->name
				}
			} sort usercmp $group->users;

			$group_rows{ $group->id } = {
				ITEM_ID      => $item_id,
				ITEM_TYPE    => $item_type,
				GROUP_ITEM   => $group->id . ':' . $conn->parent_type,
				GROUP_NAME   => $group->name,
				GROUP_ROLE   => $group->role->name,
				GROUP_DELETE => $USER->is_owner_editor( group => $group->id )
				  || $USER->is_admin,
				GROUP_USER_LOOP => \@users
			};
		}
	}

	foreach my $notebook ( values %notebooks ) {
		my %users;

		foreach my $conn ( $notebook->user_connectors ) {
			my $user = $conn->parent;
			next unless $user;
			$users{ $user->id } = {
				NOTEBOOK_USER_FULL_NAME => $user->display_name,
				NOTEBOOK_USER_NAME      => $user->name
			};
		}

		foreach my $conn ( $notebook->group_connectors ) {
			my $group = $conn->parent;
			foreach ( $group->users ) {
				$users{ $_->id } = {
					NOTEBOOK_USER_FULL_NAME => $_->display_name,
					NOTEBOOK_USER_NAME      => $_->name
				};
			}
		}

		$notebook_rows{ $notebook->id } = {
			NOTEBOOK_NAME      => $notebook->name,
			NOTEBOOK_USER_LOOP => [ values %users ]
		};
	}

	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . "Admin.tmpl" );
	$template->param(
		SHARE_DIALOG => 1,
		IS_EDITABLE  => $USER->is_admin || $isEditable,
		GROUP_LOOP =>
		  [ sort { $a->{GROUP_NAME} cmp $b->{GROUP_NAME} } values %group_rows ],
		USER_LOOP => [
			sort { $a->{USER_FULL_NAME} cmp $b->{USER_FULL_NAME} }
			  values %user_rows
		],
		NOTEBOOK_LOOP => [
			sort { $a->{NOTEBOOK_NAME} cmp $b->{NOTEBOOK_NAME} }
			  values %notebook_rows
		],
		ROLES     => get_roles('reader'),
		ITEM_ID   => $item_id,
		ITEM_TYPE => $item_type,
	);

	if ($isPublic) {
		$template->param( ACCESS_MSG => 'Everyone' );
	}

	return $template->output;
}

sub get_roles {
	my $selected = shift;

	my $html;
	foreach my $role ( $coge->resultset('Role')->all() ) {
		next
		  if ( $role->name =~ /admin/i );   # && !$USER->is_admin); # skip admin
		next if ( $role->name =~ /owner/i && !$USER->is_admin );    # skip owner
		my $name = $role->name;
		$name .= ": " . $role->description if $role->description;

#push @roles, { RID => $role->id, NAME => $name, SELECTED => ($role->id == $selected_role_id) };
		$html .=
		    '<option value="'
		  . $role->id . '" '
		  . (    $role->id eq $selected
			  || $role->name =~ /$selected/i ? 'selected' : '' )
		  . '>'
		  . $role->name
		  . '</option>';
	}
	return $html;
}

sub search_share {
	my %opts = @_;
	return if ( $USER->user_name eq 'public' );
	my $search_term = escape( $opts{search_term} );
	my $timestamp   = $opts{timestamp};

	#print STDERR "search_share $search_term $timestamp\n";

	my @results;

# Search for matching users
# $search_term = '%'.$search_term.'%';
# foreach ($coge->resultset('User')->search_literal(
# 		"user_name LIKE '$search_term' OR first_name LIKE '$search_term' OR last_name LIKE '$search_term'"))
	foreach ( $coge->resultset('User')->all ) {
		next
		  unless ( escape( $_->user_name ) =~ /$search_term/i
			|| escape( $_->display_name ) =~ /$search_term/i );
		my $label = $_->display_name . ' (' . $_->user_name . ')';
		my $value = $_->id . ':' . $node_types->{user};
		push @results, { 'label' => $label, 'value' => $value };
	}

	# Search for matching groups
	foreach ( $coge->resultset('UserGroup')->all ) {
		next unless ( escape( $_->name ) =~ /$search_term/i );
		my $label = $_->name . ' (' . $_->role->name . ' group)';
		my $value = $_->id . ':' . $node_types->{group};
		push @results, { 'label' => $label, 'value' => $value };
	}

	return encode_json( { timestamp => $timestamp, items => \@results } );
}

sub usercmp {
	no warnings 'uninitialized';    # disable warnings for undef values in sort
	$a->display_name cmp $b->display_name;
}

sub delete_genome {
	unless ( $USER->is_admin ) {
		return 0;
	}
	
	my %opts = @_;
	my $gid  = $opts{gid};
	print STDERR "delete_genome $gid\n";
	return 0 unless $gid;

	my $genome = $coge->resultset('Genome')->find($gid);
	return 0 unless $genome;
	return 0 unless ( $USER->is_admin or $USER->is_owner( dsg => $gid ) );
	my $delete_or_undelete = ( $genome->deleted ? 'undelete' : 'delete' );

	#print STDERR "delete_genome " . $genome->deleted . "\n";
	$genome->deleted( !$genome->deleted );    # do undelete if already deleted
	$genome->update;

	# Record in log
	CoGe::Accessory::Web::log_history(
		db          => $coge,
		user_id     => $USER->id,
		page        => "Admin",
		description => "$delete_or_undelete genome id $gid"
	);

	return 1;
}

sub delete_list {
	unless ( $USER->is_admin ) {
        return 0;
    }
    
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;    #return "No LID specified" unless $lid;

	my $list = $coge->resultset('List')->find($lid);
	return 0 unless $list;    #return "Cannot find list $lid\n" unless $list;
	return 0 unless ( $USER->is_admin or $USER->is_owner( list => $lid ) );
	my $delete_or_undelete = ( $list->deleted ? 'undelete' : 'delete' );

	if ( $list->locked && !$USER->is_admin ) {
		return
		  0;   #"This is a locked list.  Admin permission is needed to modify.";
	}

	$list->deleted( !$list->deleted );    # do undelete if already deleted
	$list->update;
	
	# Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => "Admin",
        description => "$delete_or_undelete notebook id $lid"
    );

	return 1;
}

sub delete_experiment {
	unless ( $USER->is_admin ) {
        return 0;
    }
    
    my %opts = @_;
    my $eid  = $opts{eid};
    return 0 unless $eid;    #return "No EID specified" unless $eid;

    my $experiment = $coge->resultset('Experiment')->find($eid);
    return 0 unless $experiment;    #return "Cannot find experiment $eid\n" unless $experiment;
    return 0 unless ( $USER->is_admin or $USER->is_owner( experiment => $eid ) );
    my $delete_or_undelete = ( $experiment->deleted ? 'undelete' : 'delete' );

    #if ( $experiment->locked && !$USER->is_admin ) {
    #    return
    #      0;   #"This is a locked experiment.  Admin permission is needed to modify.";
    #}

    $experiment->deleted( !$experiment->deleted );    # do undelete if already deleted
    $experiment->update;

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => "Admin",
        description => "$delete_or_undelete experiment id $eid"
    );

    return 1;
}

#Evan's modular stuff
#sub search_genomem {
#	# terms is an any array ference tags is a hash references
#	my ($terms, $tags, $organism_rs) = @_;
#
#	my @relevant_fields = qw(name description restricted access_count):
#
#	%organisms_tags{@relevant_fields} = $tags->{@relevant_fields};
#	$organism_tags = # ... filters tags to only organism search
#
#	my $search1 = $coge->resultset("Organism")->search({
#		-or {
#			-and {
#				create_like_fields($terms)
#			},
#			$organism_tags
#		}
#	}) ;
#
#	# Query only if organism rs exists
#	my $search2 = $organism_rs->dosomething if $organism_rs;
#
#	return join $search1 $search2;
#}

#sub search_users {
#    my %opts        = @_;
#    my $search_term = $opts{search_term};
#    my $timestamp   = $opts{timestamp};
#
#    #print STDERR "$search_term $timestamp\n";
#    return unless $search_term;
#
#    # Perform search
#    $search_term = '%' . $search_term . '%';
#    my @users = $coge->resultset("User")->search(
#        \[
#            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
#            [ 'user_name',  $search_term ],
#            [ 'first_name', $search_term ],
#            [ 'last_name',  $search_term ]
#        ]
#    );
#
#    # Limit number of results displayed
#    # if (@users > $MAX_SEARCH_RESULTS) {
#    # 	return encode_json({timestamp => $timestamp, items => undef});
#    # }
#
#    return encode_json(
#        {
#            timestamp => $timestamp,
#            items     => [ sort map { $_->user_name } @users ]
#        }
#    );
#}
