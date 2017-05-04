#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(format_time_diff);
use CoGe::JEX::Jex;
use CoGe::Core::Search;
use HTML::Template;
use JSON qw(encode_json);
use Data::Dumper;
use List::Compare;
use CoGeX;
use CoGeX::Result::User;
use CoGeDBI;
use URI::Escape::JavaScript qw(escape unescape);
use Time::Piece;
use JSON::XS;
use CoGe::Accessory::Web;
use Benchmark;
no warnings 'redefine';

$|=1;

use vars qw(
    $config $PAGE_TITLE $PAGE_NAME $user $BASEFILE $db %FUNCTION
    $FORM $MAX_SEARCH_RESULTS %ITEM_TYPE $JEX $link
);

$PAGE_TITLE = 'Admin';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $db, $user, $config, $link ) = CoGe::Accessory::Web->init( cgi => $FORM, page_title => $PAGE_TITLE );

$JEX = CoGe::JEX::Jex->new( host => $config->{JOBSERVER}, port => $config->{JOBPORT} );

$MAX_SEARCH_RESULTS = 400;

our $node_types = CoGeX::node_types();

%FUNCTION = (
	user_is_admin					=> \&user_is_admin,
	search_users                    => \&search_users,
	search_stuff                    => \&search_stuff,
	user_info                       => \&user_info,
	add_items_to_user_or_group      => \&add_items_to_user_or_group,
	remove_items_from_user_or_group => \&remove_items_from_user_or_group,
	get_share_dialog                => \&get_share_dialog,
	get_roles                       => \&get_roles,
	search_share                    => \&search_share,
	modify_item                     => \&modify_item,
    cancel_job                      => \&cancel_job,
    restart_job                     => \&restart_job,
    get_jobs                        => \&get_jobs_for_user,
    get_group_dialog                => \&get_group_dialog,
    add_users_to_group              => \&add_users_to_group,
    remove_user_from_group          => \&remove_user_from_group,
    change_group_role               => \&change_group_role,
    get_history            			=> \&get_history_for_user,
    toggle_star                     => \&toggle_star,
    update_comment                  => \&update_comment,
    update_history                  => \&update_history,
    get_user_nodes  				=> \&get_user_nodes,
    get_group_nodes  				=> \&get_group_nodes,
    get_user_table					=> \&get_user_table,
    get_group_table					=> \&get_group_table,
    get_total_table					=> \&get_total_table,
    get_jobs_table					=> \&get_jobs_table,
    get_user_jobs_table				=> \&get_user_jobs_table,
    get_user_jobs				    => \&get_user_jobs,
    gen_tree_json					=> \&gen_tree_json, 
    get_total_queries				=> \&get_total_queries,
    get_uptime		                => \&get_uptime,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
	my $html;
	my $template =
	  HTML::Template->new( filename => $config->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( USER       => $user->display_name || '',
	                  PAGE_TITLE => qq{Admin},
	                  PAGE_LINK  => $link,
	                  SUPPORT_EMAIL => $config->{SUPPORT_EMAIL},
	                  TITLE      => "GODVIEW",
	                  HOME       => $config->{SERVER},
                      HELP       => '',
                      WIKI_URL   => $config->{WIKI_URL} || '',
                      CAS_URL    => $config->{CAS_URL} || '',
                      ADMIN_ONLY => $user->is_admin,
                      COOKIE_NAME => $config->{COOKIE_NAME} || '');
	$template->param( LOGON      => 1 ) unless $user->user_name eq "public";
	$template->param( BODY       => gen_body() );
	$html .= $template->output;
}

sub gen_body {
	# Hide this page if the user is not an Admin
	unless ( $user->is_admin ) {
		my $template =
		  HTML::Template->new( filename => $config->{TMPLDIR} . "Admin.tmpl" );
		$template->param( ADMIN_ONLY => 1 );
		return $template->output;
	}

	my $template =
	  HTML::Template->new( filename => $config->{TMPLDIR} . 'Admin.tmpl' );
	$template->param( 	MAIN 			=> 1,
						API_BASE_URL  	=> $config->{SERVER} . 'api/v1/',  
						USER			=> $user->user_name,	
					);

	return $template->output;
}

sub user_is_admin {
	return $user->is_admin;
}

sub user_info {

	my %opts        = @_;
	my $search_term = $opts{search_term};
	my $search_type = $opts{search_type};
	my $timestamp   = $opts{timestamp};

	my $user;
	if ( $search_type eq "user" ) {
		$user = $db->resultset("User")->find($search_term);
	}
	if ( $search_type eq "group" ) {
		$user = $db->resultset("UserGroup")->find($search_term);
	}
	my @results;
	my @users;
	push( @users, $user );

	my $child;
	if ( $search_type eq "user" ) {
		foreach ( $user->child_connectors( { child_type => 6 } ) ) {
			$child = $_->child;
			push( @users, $child );
		}
	}

	my $i;
	for ($i = 0; $i < scalar(@users); $i++) {
		my $currentUser = $users[$i];
		my @current;

		# Find notebooks
		foreach ( $currentUser->child_connectors( { child_type => 1 } ) ) {
			$child = $_->child;
			if ($user->has_access_to_notebook($child)) {
				push @current,
				  {
					'type'          => "notebook",
					'label'         => $child->info,
					'id'            => $child->id,
					'role'          => $_->role_id,
					'deleted'       => $child->deleted,
					'restricted'    => $child->restricted
				  };
			}
		}

		# Find genomes
		foreach ( $currentUser->child_connectors( { child_type => 2 } ) ) {
			$child = $_->child;
			if ($user->has_access_to_genome($child)) {
				push @current,
				  {
					'type'          => "genome",
					'label'         => $child->info,
					'id'            => $child->id,
					'role'          => $_->role_id,
					'deleted'       => $child->deleted,
	                'restricted'    => $child->restricted
				  };
			}
		}

		# Find experiments
		foreach ( $currentUser->child_connectors( { child_type => 3 } ) ) {
			$child = $_->child;
			if ($user->has_access_to_experiment($child)) {
				push @current,
				  {
					'type'          => "experiment",
					'label'         => $child->name,
					'id'            => $child->id,
					'info'          => $child->info,
					'role'          => $_->role_id,
					'deleted'       => $child->deleted,
	                'restricted'    => $child->restricted
				  };
			}
		}

		# Find users if searching a user group
		if ( $search_type eq "group" ) {	
			foreach ( $user->users ) {
				push @current,
				  { 'type' => "user", 'first' => $_->first_name, 'last' => $_->last_name, 'username' => $_->user_name, 'id' => $_->id, 'email' => $_->email };
			}
		}
		
		if ($search_type eq "user" && $i == 0) {
			push @results,
			  {
			  	'first'		=> $currentUser->first_name,
			  	'last'		=> $currentUser->last_name,
				'username'	=> $currentUser->user_name,
				'email'		=> $currentUser->email,
				'user_id' 	=> $currentUser->id,
				'result' 	=> \@current
			  };	  
		} else {
			push @results,
			  {
				'user'		=> $currentUser->name,
				'user_id' 	=> $currentUser->id,
				'result' 	=> \@current
			  };
		}
	}

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

	my ( $target_id, $target_type ) = $target_item =~ /(\d+)\:(\d+)/;

	return unless ( $target_id and $target_type );

	# Verify that user has access to each item
	my @verified;
	foreach my $item (@items) {
		my ( $item_id, $item_type ) = $item =~ /content_(\d+)_(\d+)/;

		next unless ( $item_id and $item_type );

		if ( $item_type == $node_types->{genome} ) {
			my $genome = $db->resultset('Genome')->find($item_id);
			next
			  unless ( $user->has_access_to_genome($genome)
				|| $user->is_admin );
			push @verified, $item;
		}
		elsif ( $item_type == $node_types->{experiment} ) {
			my $experiment = $db->resultset('Experiment')->find($item_id);
			next unless $user->has_access_to_experiment($experiment);
			push @verified, $item;
		}
		elsif ( $item_type == $node_types->{notebook} ) {
			my $notebook = $db->resultset('List')->find($item_id);
			next unless $user->has_access_to_notebook($notebook);
			push @verified, $item;
		}
	}

	# Assign each item to user/group
	#TODO verify that user can use specified role (for admin/owner roles)
	if ( $target_type == $node_types->{user} ) {
		my $user = $db->resultset('User')->find($target_id);

		return unless $user;

		foreach (@verified) {
			my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

			# Remove previous connection
			foreach (
				$db->resultset('UserConnector')->search(
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
			my $conn = $db->resultset('UserConnector')->create(
				{
					parent_id   => $target_id,
					parent_type => 5,            # FIXME hardcoded
					child_id    => $item_id,
					child_type  => $item_type,
					role_id     => $role_id
				}
			);

			return unless $conn;
		}
	}
	elsif ( $target_type == $node_types->{group} ) {
		my $group = $db->resultset('UserGroup')->find($target_id);

		return unless $group;

		foreach (@verified) {
			my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

			my $conn = $db->resultset('UserConnector')->find_or_create(
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

	return unless $target_item;
	my $item_list = $opts{item_list};
	my @items = split( ',', $item_list );

	return unless @items;

	my ( $target_id, $target_type ) = $target_item =~ /(\d+)\:(\d+)/;

	next unless ( $target_id and $target_type );

	foreach (@items) {
		my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

		next unless ( $item_id and $item_type );

		if ( $item_type == $node_types->{genome} ) {
			my $genome = $db->resultset('Genome')->find($item_id);
			next unless ( $user->has_access_to_genome($genome) );

			my $conn = $db->resultset('UserConnector')->find(
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
			my $experiment = $db->resultset('Experiment')->find($item_id);
			next unless $user->has_access_to_experiment($experiment);

			my $conn = $db->resultset('UserConnector')->find(
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
			my $notebook = $db->resultset('List')->find($item_id);
			next unless $user->has_access_to_notebook($notebook);

			my $conn = $db->resultset('UserConnector')->find(
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
	my %opts      = @_;
	my $item_list = $opts{item_list};

	my @items = split( ',', $item_list );
	return unless @items;

	my ( %userconn, %notebooks );
	my $isPublic   = 0;
	my $isEditable = 1;
	my $item_id;
	my $item_type;
	foreach (@items) {
		( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

		next unless ( $item_id and $item_type );

		if ( $item_type == $node_types->{genome} ) {
			my $genome = $db->resultset('Genome')->find($item_id);
			next unless ( $user->has_access_to_genome($genome) );
			map { $userconn{ $_->parent_id } = $_ }
			  ( $genome->user_connectors, $genome->group_connectors );
			map { $notebooks{ $_->id } = $_ } $genome->lists;
			$isPublic = 1 if ( not $genome->restricted );
			$isEditable = 0 if ( not $user->is_owner_editor( dsg => $genome ) );
		}
		elsif ( $item_type == $node_types->{experiment} ) {
			my $experiment = $db->resultset('Experiment')->find($item_id);
			next unless $user->has_access_to_experiment($experiment);
			map { $userconn{ $_->id } = $_ }
			  ( $experiment->user_connectors, $experiment->group_connectors );
			map { $notebooks{ $_->id } = $_ } $experiment->lists;
			$isPublic = 1 if ( not $experiment->restricted );
			$isEditable = 0
			  if ( not $user->is_owner_editor( experiment => $experiment ) );
		}
		elsif ( $item_type == $node_types->{notebook} ) {
			my $notebook = $db->resultset('List')->find($item_id);
			next unless $user->has_access_to_notebook($notebook);
			map { $userconn{ $_->id } = $_ }
			  ( $notebook->user_connectors, $notebook->group_connectors );
			$isPublic = 1 if ( not $notebook->restricted );
			$isEditable = 0
			  if ( not $user->is_owner_editor( list => $notebook ) );
		}
	}

	$isEditable = 1 if ( $user->is_admin );

	my ( %user_rows, %group_rows, %notebook_rows );
	foreach my $conn ( values %userconn ) {
		next
		  unless $conn->role
		; #EL added 10/21/14 so solve a problem for a user where the share dialog wouldn't appear for his genomes, and an fatal error was being thrown due to role->name with role being undefined.  Genomes with problems were IDs: 24576, 24721, 24518, 24515, 24564, 24566, 24568, 24562, 24571
		if ( $conn->is_parent_user ) {
			my $user = $conn->parent;

			$user_rows{ $user->id } = {
				ITEM_ID        => $item_id,
				ITEM_TYPE      => $item_type,
				USER_ITEM      => $user->id . ':' . $conn->parent_type,
				USER_FULL_NAME => $user->display_name,
				USER_NAME      => $user->name,
				USER_ROLE      => $conn->role->name,
				USER_DELETE    => $isEditable
				  && (!$conn->role->is_owner
					|| $user->is_admin )   # owner can't be removed unless admin
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
				GROUP_DELETE => $user->is_owner_editor( group => $group->id )
				  || $user->is_admin,
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
	  HTML::Template->new( filename => $config->{TMPLDIR} . "Admin.tmpl" );
	$template->param(
		SHARE_DIALOG => 1,
		IS_EDITABLE  => $user->is_admin || $isEditable,
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

sub get_group_dialog {
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;
    my $item_id;
    my $item_type;

    my ( %users, %roles, %creators, %owners, $lowest_role );
    foreach (@items) {
        ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
        next unless ( $item_type == $node_types->{group} ); # sanity check

        my $group = $db->resultset('UserGroup')->find($item_id);
        next unless ( $group and $group->is_editable($user) );
        next if ( $group->locked and !$user->is_admin );

        my $role = $group->role;
        $lowest_role = $role if (!$lowest_role or $role->is_lower($lowest_role));
        my $creator_id = $group->creator_user_id;
        my $owner_id = $group->owner->id;
        foreach my $user ($group->users) {
            my $uid = $user->id;
            $users{$uid} = $user;
            if ( not defined $roles{$uid}
                 or $role->is_lower($roles{$uid}) )
            {
                $roles{$uid} = $role;
            }
            if ($uid == $creator_id) {
                $creators{$uid} = 1;
            }
            if ($uid == $owner_id) {
                $owners{$uid} = 1;
            }
        }
    }

    my @rows;
    foreach my $user ( sort usercmp values %users ) {
        my $uid = $user->id;
        my $role_name;
        if ($creators{$uid}) {
            $role_name = 'Creator';
        }
        if ($owners{$uid}) {
            $role_name = ($role_name ? $role_name . ', ' : '') . 'Owner';
        }
        push @rows, {
        	ITEM_ID      => $item_id,
            ITEM_TYPE    => $item_type,
            USER_ITEM      => $uid,
            USER_FULL_NAME => $user->display_name,
            USER_NAME      => $user->name,
            USER_ROLE      => ($role_name ? ' - ' . $role_name : ''),#$roles{$uid}->name),
            USER_DELETE    => !$owners{$uid} # owner can't be removed
        };
    }

    my $template =
      HTML::Template->new( filename => $config->{TMPLDIR} . "Admin.tmpl" );

    # If no editable groups then show error dialog
    if (!$user->is_admin) {
        $template->param(
            ERROR_DIALOG => 1,
            ERROR_MESSAGE => "You don't have permission to modify the selected group(s).",
        );
    }
    else {
        $template->param(
            GROUP_DIALOG => 1,
            IS_EDITABLE  => 1,
            USER_LOOP => \@rows,
            ROLES => get_roles($lowest_role->name),
            ITEM_ID   => $item_id,
            ITEM_TYPE => $item_type,
        );
    }

    return $template->output;
}

sub get_roles {
	my $selected = shift;

	my $html;
	foreach my $role ( $db->resultset('Role')->all() ) {
		next
		  if ( $role->name =~ /admin/i );   # && !$user->is_admin); # skip admin
		next if ( $role->name =~ /owner/i && !$user->is_admin );    # skip owner
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
	return if ( $user->user_name eq 'public' );
	my $search_term = escape( $opts{search_term} );
	my $timestamp   = $opts{timestamp};

	my @results;

# Search for matching users
# $search_term = '%'.$search_term.'%';
# foreach ($db->resultset('User')->search_literal(
# 		"user_name LIKE '$search_term' OR first_name LIKE '$search_term' OR last_name LIKE '$search_term'"))
	foreach ( $db->resultset('User')->all ) {
		next
		  unless ( escape( $_->user_name ) =~ /$search_term/i
			|| escape( $_->display_name ) =~ /$search_term/i );
		my $label = $_->display_name . ' (' . $_->user_name . ')';
		my $value = $_->id . ':' . $node_types->{user};
		push @results, { 'label' => $label, 'value' => $value };
	}

	# Search for matching groups
	foreach ( $db->resultset('UserGroup')->all ) {
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

sub modify_item {
	
	unless ( $user->is_admin ) {
        return 0;
    }
    
    my %opts = @_;
    my $id  = $opts{id};
    my $mod  = $opts{modification};
    my $type = $opts{type};
    return 0 unless $id;
    my $item_type = $node_types->{lc($type)}; # mdb added 7/21/15 for newsfeed
    die unless $item_type;

    my $item = $db->resultset($type)->find($id);
    return 0 unless $item;
    return 0 unless ( $user->is_admin or $user->is_owner( dsg => $id ) );
    
    my $log_message;
    if ($mod eq "delete") {
        $log_message = ( $item->deleted ? 'undeleted' : 'deleted' );
        $item->deleted( !$item->deleted );    # do undelete if already deleted
    } 
    elsif ($mod eq "restrict") {
    	$log_message = ( $item->restricted ? 'unrestricted' : 'restricted' );
    	$item->restricted( !$item->restricted );
    }
    $item->update;

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $db,
        user_id     => $user->id,
        page        => "Admin",
        description => "$log_message $type " . $item->info_html,
        parent_id   => $id,
        parent_type => $item_type
    );

    return 1;
}

sub add_users_to_group {
    my %opts = @_;
    my @target_items = split( ',', $opts{target_items} );
    return unless @target_items;
    my $new_item = $opts{new_item};
    return unless $new_item;

    # Build a list of users to add to the target group
    my %users;
    my ( $item_id, $item_type ) = $new_item =~ /(\d+)\:(\d+)/;
    return unless ( $item_id and $item_type );

    if ( $item_type == $node_types->{user} ) {
        my $user = $db->resultset('User')->find($item_id);
        return unless $user;
        $users{$user->id} = $user;
    }
    elsif ( $item_type == $node_types->{group} ) {
        my $group = $db->resultset('UserGroup')->find($item_id);
        return unless $group;
        # TODO check that user has visibility of this group (one that they own or belong to)
        map { $users{$_->id} = $_ } $group->users;
     }

    # Add users to the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
        my ( $target_id, $target_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $target_id and $target_type );
        next unless ( $target_type == $node_types->{group} ); # sanity check
        my $target_group = $db->resultset('UserGroup')->find($target_id);
        next unless ( $target_group and $target_group->is_editable($user) );
        next if ( $target_group->locked && !$user->is_admin );

        # Add users to this target group
        foreach my $user (values %users) {
            # Check for existing user connection to target group
            my $conn = $db->resultset('UserConnector')->find(
                {
                    parent_id   => $user->id,
                    parent_type => 5,                #FIXME hardcoded to "user"
                    child_id    => $target_id,
                    child_type  => 6                 #FIXME hardcoded to "group"
                }
            );

            # Create new user connection if one wasn't found
            if (!$conn) {
                $conn = $db->resultset('UserConnector')->create(
                    {
                        parent_id   => $user->id,
                        parent_type => 5,                #FIXME hardcoded to "user"
                        child_id    => $target_id,
                        child_type  => 6,                #FIXME hardcoded to "group"
                        role_id     => $target_group->role_id
                    }
                );
            }
            next unless $conn;

            # Record in log
            CoGe::Accessory::Web::log_history(
                db          => $db,
                user_id     => $user->id,
                page        => "Admin",
                description => 'added user ' . $user->info . ' to group ' . $target_group->info_html,
                parent_id   => $target_id,
                parent_type => 6 #FIXME magic number
            );
        }
    }

    return get_group_dialog( item_list => $opts{target_items} );
}

sub remove_user_from_group {
    my %opts = @_;
    my @target_items = split( ',', $opts{target_items} );
    return unless @target_items;
    my $user_id = $opts{user_id};
    return unless $user_id;

    # Verify user
    my $user = $db->resultset('User')->find($user_id);
    return unless $user;

    # Remove users from the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
        my ( $target_id, $target_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $target_id and $target_type );
        next unless ( $target_type == $node_types->{group} ); # sanity check
        my $target_group = $db->resultset('UserGroup')->find($target_id);
        next unless ( $target_group and $target_group->is_editable($user) );
        next if ( $target_group->locked && !$user->is_admin );

        # Get user connection to target group
        my $conn = $db->resultset('UserConnector')->find(
            {
                parent_id   => $user_id,
                parent_type => 5,                #FIXME hardcoded to "user"
                child_id    => $target_id,
                child_type  => 6                 #FIXME hardcoded to "group"
            }
        );
        next unless $conn;

        # Delete user connection if not owner
        next if ($conn->role->is_owner);
        $conn->delete;

        # Record in log
        CoGe::Accessory::Web::log_history(
            db          => $db,
            user_id     => $user->id,
            page        => "Admin",
            description => 'removed user ' . $user->info . ' from group ' . $target_group->info_html,
            parent_id   => $target_id,
            parent_type => 6
        );
    }

    return get_group_dialog( item_list => $opts{target_items} );
}

sub change_group_role {
    my %opts = @_;
    my @target_items = split( ',', $opts{target_items} );
    return unless @target_items;
    my $role_id = $opts{role_id};
    return unless $role_id;

    # Verify role
    my $role = $db->resultset('Role')->find($role_id);
    return unless $role;

    # Change role for the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
        my ( $target_id, $target_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $target_id and $target_type );
        next unless ( $target_type == $node_types->{group} ); # sanity check
        my $target_group = $db->resultset('UserGroup')->find($target_id);
        next unless ( $target_group and $target_group->is_editable($user) );
        next if ( $target_group->locked && !$user->is_admin );

        $target_group->role_id($role_id);
        $target_group->update;
    }

    return get_group_dialog( item_list => $opts{target_items} );
}

sub get_jobs_for_user {
    my %opts = @_;
	my $status 		 = $opts{status}; # running, completed, failed, etc...
	my $name         = $opts{name};   # SynMap, CoGeBlast, LoadGenome, etc...

	return encode_json({ error => 'Not logged in' }) if ($user->is_public);

	# Get workflows from DB
	my $where = 'parent_id IS NOT NULL and type != 0';
	$where .= ' AND page=' . $db->dbh->quote($name) if $name;
	$where .= ' AND log.user_id=' . $user->id if !$user->is_admin;
	my $logs = $db->storage->dbh->selectall_arrayref('SELECT page,link,parent_id,user_name FROM (SELECT page,link,parent_id,user_name,time FROM log JOIN user ON user.user_id=log.user_id WHERE ' . $where . ' ORDER BY log_id DESC) t GROUP BY parent_id ORDER BY time DESC');

	# Get workflows from JEX
    my $workflows;
    if ($status) { # Get workflows with given status
		$workflows = $JEX->find_workflows(undef, $status);
	}
	else { # Get only workflows from DB query
		my @workflow_ids = map { $_->[2] } @{$logs};
    	$workflows = $JEX->find_workflows(\@workflow_ids);
    }

    my %workflow_results;
    foreach (@{$workflows}) {
        my($id, $name, $submitted, $completed, $status) = @{$_};

        my $start_time = localtime($submitted)->strftime('%F %H:%M');
        my $end_time = "";
        my $diff;

        if ($completed) {
            $end_time = localtime($completed)->strftime('%F %H:%M');
            $diff = $completed - $submitted;
        }
		else {
            $diff = time - $submitted;
        }

        $workflow_results{$id} = {
            status    => $status,
            started   => $start_time,
            completed => $end_time,
            elapsed   => format_time_diff($diff)
        };
    }

	# Intersect DB log with JEX workflows
	my @job_items;
    foreach (@{$logs}) {
        my $entry = $workflow_results{int($_->[2])};
        next unless $entry; # a log entry must correspond to a workflow
		# [ID, Started, Completed, Elapsed, User, Tool, Link, Status]
        push @job_items, [
            $_->[2], # parent_id
            $entry->{started},
            $entry->{completed},
            $entry->{elapsed},
            $_->[3] || "public",
            $_->[0], # page
            $_->[1], # link
            $entry->{status},
        ];
    }

    # Remove repeated entries
	# my @filtered;
    # foreach (reverse @job_items) {
    # 	my @job = @$_;
    #     my $wid = $job[0];
    #     next if (defined $wid and defined $workflow_results{$wid}{seen});
    #     $workflow_results{$wid}{seen}++ if (defined $wid);
    #     unshift @filtered, \@job;
    # }

    return encode_json({ # formatted for DataTables
    	data => \@job_items,
    	#bPaginate => 0,
		columnDefs => [{ 
			orderSequence => [ "desc", "asc" ], 
			targets => [0, 1, 2],
		}],
    });
}

sub cancel_job {
    my $job_id = _check_job_args(@_);

    return encode_json( {} ) unless defined($job_id);

    my $status = $JEX->get_status( $job_id );

    if ( $status =~ /scheduled|running|notfound/i ) {
        return encode_json({ status => $JEX->terminate( $job_id ) });
    } else {
        return encode_json( {} );
    }
}

sub restart_job {
    my $job_id = _check_job_args(@_);

    return encode_json( {} ) unless defined($job_id);

    my $status = $JEX->get_status( $job_id );

    if ( $status =~ /running/i ) {
        return encode_json( {} );
    } else {
        return encode_json({ status => $JEX->restart( $job_id ) });
    }
}

sub cmp_by_start_time {
    my $job1 = shift;
    my $job2 = shift;

    $job1->start_time cmp $job2->start_time;
}

sub _check_job_args {
    my %args   = @_;
    my $job_id = $args{job};

    if ( not defined($job_id) ) {
        say STDERR "Job.pl: a job id was not given to cancel_job.";
    }

    return $job_id;
}


#History Tab
sub get_history_for_user {
    my %opts       = @_;
    my $time_range = $opts{time_range};    # in hours
    my $timestamp = $opts{timestamp};
    
    if ($timestamp) {
	    $timestamp = "'" . $timestamp . "'";
    } else {
    	$timestamp = "'9999-12-31 23:59:59'";
    }
    
    $time_range = 24 if ( not defined $time_range or $time_range !~ /[-\d]/ );
    my $include_page_accesses = $opts{include_pages_accesses};

    my %users = map { $_->user_id => $_->name } $db->resultset('User')->all;
    
    my @entries;
	if ( $user->is_admin ) {
		if ( $time_range == 0 ) {
			@entries = $db->resultset('Log')->search(
				{ 
					'time' => { '<' => \$timestamp }, 
					type => { '!=' => 0 }
				},
				{ 
					order_by => { -desc => 'time' },
					rows => 1000
				}
			);
		}
	}
	else {
		if ( $time_range == 0 or $time_range == -3 ) {
			@entries = $db->resultset('Log')->search(
				{
					user_id => $user->id,
					'time' => { '<' => \$timestamp },
					#{ description => { 'not like' => 'page access' } }
					type => { '!=' => 0 }
				},
				{ 
					order_by => { -desc => 'time' },
					rows => 1000
				}
			);
		}
	}

    my @items;
    foreach (@entries) {
	    #[Date/Time, User, Page, Descrition, Link, Comment]
	    push @items, [
			#id          => $_->id,
			#starred     => ( $_->status != 0 ),
			$_->time,
			( $_->user_id ? $users{ $_->user_id } : 'public' ),
			$_->page,
			$_->description,
			( $_->link ? $_->link : '' ),
			$_->comment
		];
    }

    return encode_json({data => \@items});
}

sub toggle_star {
    my %opts   = @_;
    my $log_id = $opts{log_id};

    my $entry = $db->resultset('Log')->find($log_id);
    return '' unless $entry;

    my $status = $entry->status;
    $entry->status( not $status );
    $entry->update();

    return not $status;
}

sub update_comment {
    my %opts    = @_;
    my $log_id  = $opts{log_id};
    my $comment = $opts{comment};

    my $entry = $db->resultset('Log')->find($log_id);
    return unless $entry;

    $entry->comment($comment);
    $entry->update();
}

sub update_history {
	my %opts       = @_;
    my $time_range = $opts{time_range};    # in hours
    my $timestamp = $opts{timestamp};
    $timestamp = "'" . $timestamp . "'";
    
    $time_range = 24 if ( not defined $time_range or $time_range !~ /[-\d]/ );
    my $include_page_accesses = $opts{include_pages_accesses};

    my %users = map { $_->user_id => $_->name } $db->resultset('User')->all;
    my @entries;
    if ( $user->is_admin ) {
        if ( $time_range == 0 ) {
            @entries = $db->resultset('Log')->search(
                { 'time' => { '>' => \$timestamp } },
                { order_by => { -desc => 'time' } }
            );
        }
    }
    
    my @items;
    foreach (@entries) {
    	#[Date/Time, User, Page, Descrition, Link, Comment]
        push @items, [
            #id          => $_->id,
            #starred     => ( $_->status != 0 ),
			$_->time,
			( $_->user_id ? $users{ $_->user_id } : 'public' ),
			$_->page,
			$_->description,
			( $_->link ? $_->link : '' ),
			$_->comment
		];
    }

    return encode_json({new_rows => \@items});
}	


#Graph tab
sub get_user_nodes {
    my %childrenByList;
    foreach my $conn ( $db->resultset('ListConnector')->all ) {
    	if($conn->child_type != 4) {
    		my $child = $conn->child;
    		
    		if($child) {
	        	push @{ $childrenByList{ $conn->parent_id } },
          		{ 
          			name 		=> $conn->child_id, 
          			size 		=> 2025, 
          			type 		=> $conn->child_type,
          			deleted		=> $child->deleted,
          			restricted 	=> $child->restricted,
          		};
    		}
    	}
    }

    my %childrenByUser;
    foreach my $conn (
        $db->resultset('UserConnector')->search(
            {
                parent_type => $node_types->{user},
                -or         => [
                    child_type => $node_types->{genome},
                    child_type => $node_types->{experiment},
                    child_type => $node_types->{list}
                ]
            }
        )
      )
    {
        if ( $conn->child_type == $node_types->{list} ) {
        	my @children = $childrenByList{ $conn->child_id};
    		my $sub = $children[0];
    		my $size = 1;
    		if($sub) {
    			$size = scalar @{$sub};
    		}
    		
    		my $child = $conn->child;
    		
            push @{ $childrenByUser{ $conn->parent_id } },
            {
           		name     	=> $conn->child_id,
                type     	=> $conn->child_type,
                size	 	=> $size * 1000,
                _size	 	=> 2025,
                children 	=> $childrenByList{ $conn->child_id },
                deleted	 	=> $child->deleted,
                restricted 	=> $child->restricted,
              };
        }
        else {
        	my $child = $conn->child;
            push @{ $childrenByUser{ $conn->parent_id } },
              { 
              	name 		=> $conn->child_id, 
              	size 		=> 2025, 
              	type 		=> $conn->child_type,
              	deleted 	=> $child->deleted,
              	restricted 	=> $child->restricted,
              };
        }
    }

    my @users;
    foreach my $user ( $db->resultset('User')->all ) {
    	my @children = $childrenByUser{ $user->id};
    	my $sub = $children[0];
    	my $size = 1;
    	if($sub) {
    		$size = scalar @{$sub};
    	}
    	
        push @users,
          {
            name     	=> $user->id,
            type     	=> 5,
            info     	=> $user->info,
            size	 	=> $size * 1000,
            _size	 	=> 2025,
            children 	=> $childrenByUser{ $user->id },
            deleted	 	=> 0,
            restricted 	=> 0,
          };
    }    

    return encode_json(
		{ 
			name 	 => 'users',  
			info	 => 'users', 
			children => \@users, 
			size	 => (scalar @users) * 1000, 
			_size	 => 2025,
		}
    );
}

sub get_group_nodes {
	my %childrenByList;
	# TODO: Only pull from relevant lists
    foreach my $conn ( $db->resultset('ListConnector')->all ) {
    	if($conn->child_type != 4) {
    		my $child = $conn->child;
    		
    		if($child) {
	        	push @{ $childrenByList{ $conn->parent_id } },
          		{ 
          			name 		=> $conn->child_id, 
          			size 		=> 2025, 
          			type 		=> $conn->child_type,
          			deleted		=> $child->deleted,
          			restricted 	=> $child->restricted,
          		};
    		}
    	}
    }
	
	my %childrenByGroup;
    foreach my $conn (
        $db->resultset('UserConnector')->search(
            {
                parent_type => $node_types->{group},
                -or         => [
                    child_type => $node_types->{genome},
                    child_type => $node_types->{experiment},
                    child_type => $node_types->{list}
                ]
            }
        )
      )
    {
        if ( $conn->child_type == $node_types->{list} ) {
        	my @children = $childrenByList{ $conn->child_id};
    		my $sub = $children[0];
    		my $size = 1;
    		if($sub) {
    			$size = scalar @{$sub};
    		}
    		
    		my $child = $conn->child;
    		
            push @{ $childrenByGroup{ $conn->parent_id } },
              {
                name     	=> $conn->child_id,
                type     	=> $conn->child_type,
                size	 	=> $size * 1000,
                _size	 	=> 2025,
                children 	=> $childrenByList{ $conn->child_id },
                deleted	 	=> $child->deleted,
                restricted 	=> $child->restricted,
              };
        }
        else {
        	my $child = $conn->child;
            push @{ $childrenByGroup{ $conn->parent_id } },
              { name 		=> $conn->child_id, 
              	size 		=> 2025, 
              	type 		=> $conn->child_type,
              	deleted 	=> $child->deleted,
              	restricted 	=> $child->restricted,
              };
        }
    }

    my @groups;
    foreach my $group ( $db->resultset('UserGroup')->all ) {
    	my @children = $childrenByGroup{ $group->id};
    	my $sub = $children[0];
    	my $size = 1;
    	if($sub) {
    		$size = scalar @{$sub};
    	}
        #my $num_users = @{ $group->users };
        push @groups,
          {
            name     	=> $group->id,
            type     	=> 6,
            info     	=> $group->info,
            size   	 	=> $size * 3000,
            _size	 	=> 2025,
            children 	=> $childrenByGroup{ $group->id },
            deleted	 	=> 0,
            restricted 	=> 0,
          };
    }
    
    return encode_json(
		{ 
			name	 => 'groups', 
			info	 => 'groups', 
			children => \@groups, 
			size	 => (scalar @groups) * 3000,
			_size	 => 2025,
		}
    );
}

sub _get_counts_by_user {
	my $child_type = shift;
	my $join = shift;
	my $where = shift;
	my $select = 'SELECT parent_id,FORMAT(COUNT(*), 0) FROM user_connector JOIN ' . $join . ' ON ' . $join . '_id=child_id WHERE parent_type=5 AND child_type=' . $child_type . ' AND ' . $where . ' GROUP BY parent_id';
	my $rows = $db->storage->dbh->selectall_arrayref($select);
	my %map = map { $_->[0] => $_->[1] } @$rows;
	return \%map;
}
	
sub get_user_table {
	my %opts = @_;
    my $filter = $opts{filter};
    my @data;

	my $user_rows = $db->storage->dbh->selectall_arrayref('SELECT user_id,first_name,last_name,user_name FROM user ORDER BY first_name,last_name');
	my ($g, $e, $n, $ug);

    if ($filter eq "restricted") {
		$g = _get_counts_by_user 2, 'genome', 'user_connector.role_id=2 AND deleted=0 AND restricted=1';
		$e = _get_counts_by_user 3, 'experiment', 'user_connector.role_id=2 AND deleted=0 AND restricted=1';
		$n = _get_counts_by_user 1, 'list', 'user_connector.role_id=2 AND deleted=0 AND restricted=1';
    } elsif ($filter eq "deleted") {
		$g = _get_counts_by_user 2, 'genome', 'user_connector.role_id=2 AND deleted=1';
		$e = _get_counts_by_user 3, 'experiment', 'user_connector.role_id=2 AND deleted=1';
		$n = _get_counts_by_user 1, 'list', 'user_connector.role_id=2 AND deleted=1';
		$ug = _get_counts_by_user 6, 'user_group', 'deleted=1';
    } elsif ($filter eq "public") {
		$g = _get_counts_by_user 2, 'genome', 'user_connector.role_id=2 AND deleted=0 AND restricted=0';
		$e = _get_counts_by_user 3, 'experiment', 'user_connector.role_id=2 AND deleted=0 AND restricted=0';
		$n = _get_counts_by_user 1, 'list', 'user_connector.role_id=2 AND deleted=0 AND restricted=0';
		$ug = _get_counts_by_user 6, 'user_group', 'deleted=0';
    } elsif ($filter eq "shared") {
		$g = _get_counts_by_user 2, 'genome', 'deleted=0 AND restricted=1 AND user_connector.role_id!=2';
		$e = _get_counts_by_user 3, 'experiment', 'deleted=0 AND restricted=1 AND user_connector.role_id!=2';
		$n = _get_counts_by_user 1, 'list', 'deleted=0 AND restricted=1 AND user_connector.role_id!=2';
    } else {
		$g = _get_counts_by_user 2, 'genome', 'user_connector.role_id=2 AND deleted=0';
		$e = _get_counts_by_user 3, 'experiment', 'user_connector.role_id=2 AND deleted=0';
		$n = _get_counts_by_user 1, 'list', 'user_connector.role_id=2 AND deleted=0';
		$ug = _get_counts_by_user 6, 'user_group', 'deleted=0';
	}

	my @user_data;
	foreach my $u (@$user_rows) {
		my $user_id = $u->[0];
		my $num_genomes = $g->{$user_id};
		my $num_experiments = $e->{$user_id};
		my $num_notebooks = $n->{$user_id};
		my $num_groups = $ug ? $ug->{$user_id} : undef;
		if ($num_genomes + $num_experiments + $num_notebooks > 0) {
			push @data, [$u->[1] . ' ' . $u->[2] . ' (' . $u->[3] . ': <a href="User.pl?user_id=' . $user_id . '" target="_blank">' . $user_id . '</a>)', $num_genomes, $num_experiments, $num_notebooks, $num_groups];
		}
	}
	return encode_json({
		data => \@data,
		bPaginate => 0,
		columnDefs => [{ 
			orderSequence => [ "desc", "asc" ], 
			targets => [1, 2, 3, 4, 5],
		}]
	});
}

sub _get_counts_by_group {
	my $child_type = shift;
	my $join = shift;
	my $where = shift;
	my $select = 'SELECT parent_id,FORMAT(COUNT(*), 0) FROM user_connector JOIN ' . $join . ' ON ' . $join . '_id=child_id WHERE parent_type=6 AND child_type=' . $child_type . ' AND ' . $where . ' GROUP BY parent_id';
	my $rows = $db->storage->dbh->selectall_arrayref($select);
	my %map = map { $_->[0] => $_->[1] } @$rows;
	return \%map;
}
	
sub get_group_table {
	my %opts = @_;
    my $filter = $opts{filter};
	my @data;

	my $group_rows = $db->storage->dbh->selectall_arrayref('SELECT user_group_id,name FROM user_group ORDER BY name');
	my ($g, $e, $n);

    if ($filter eq "restricted") {
		$g = _get_counts_by_group 2, 'genome', 'deleted=0 AND restricted=1';
		$e = _get_counts_by_group 3, 'experiment', 'deleted=0 AND restricted=1';
		$n = _get_counts_by_group 1, 'list', 'deleted=0 AND restricted=1';
    } elsif ($filter eq "deleted") {
		$g = _get_counts_by_group 2, 'genome', 'deleted=1';
		$e = _get_counts_by_group 3, 'experiment', 'deleted=1';
		$n = _get_counts_by_group 1, 'list', 'deleted=1';
    } elsif ($filter eq "public") {
		$g = _get_counts_by_group 2, 'genome', 'deleted=0 AND restricted=0';
		$e = _get_counts_by_group 3, 'experiment', 'deleted=0 AND restricted=0';
		$n = _get_counts_by_group 1, 'list', 'deleted=0 AND restricted=0';
    } elsif ($filter eq "shared") {
		$g = _get_counts_by_group 2, 'genome', 'deleted=0 AND restricted=1';
		$e = _get_counts_by_group 3, 'experiment', 'deleted=0 AND restricted=1';
		$n = _get_counts_by_group 1, 'list', 'deleted=0 AND restricted=1';
    } else {
		$g = _get_counts_by_group 2, 'genome', 'deleted=0';
		$e = _get_counts_by_group 3, 'experiment', 'deleted=0';
		$n = _get_counts_by_group 1, 'list', 'deleted=0';
	}

	my @group_data;
	foreach my $group (@$group_rows) {
		my $group_id = $group->[0];
		my $num_genomes = $g->{$group_id};
		my $num_experiments = $e->{$group_id};
		my $num_notebooks = $n->{$group_id};
		my $num_users = $db->storage->dbh->selectall_arrayref('SELECT FORMAT(COUNT(*), 0) FROM user_connector WHERE child_id=' . $group_id . ' AND parent_type=5 and child_type=6')->[0][0];
		if ($num_genomes + $num_experiments + $num_notebooks > 0) {
			push @data, [$group->[1] . ' (<a href="#" onclick="group_dialog(' . $group_id . ')">' . $group_id . '</a>)', $num_genomes, $num_experiments, $num_notebooks, $num_users];
		}
	}
	return encode_json({
		data => \@data,
		bPaginate => 0,
		columnDefs => [{ 
			orderSequence => [ "desc", "asc" ], 
			targets => [1, 2, 3, 4, 5],
		}]
	});
}

sub get_total_table {
	my @data;

	my $total_genomes = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM genome WHERE deleted=0');
	my $total_experiments = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM experiment WHERE deleted=0');
	my $total_notebooks = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM list WHERE deleted=0');
	my $total_users = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM user');
	my $total_groups = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM user_group WHERE deleted=0');
	push @data, [undef, $total_genomes, $total_experiments, $total_notebooks, $total_users, $total_groups];

	my $restricted_genomes = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM genome WHERE restricted=1');
	my $restricted_experiments = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM experiment WHERE restricted=1');
	my $restricted_notebooks = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM list WHERE restricted=1');
	push @data, ['restricted', $restricted_genomes, $restricted_experiments, $restricted_notebooks, undef, undef];
	
	my $deleted_genomes = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM genome WHERE deleted=1');
	my $deleted_experiments = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM experiment WHERE deleted=1');
	my $deleted_notebooks = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM list WHERE deleted=1');
	my $deleted_groups = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM user_group WHERE deleted=1');
	push @data, ['deleted', $deleted_genomes, $deleted_experiments, $deleted_notebooks, undef, $deleted_groups];
	
	my $public_genomes = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM genome WHERE restricted=0 AND deleted=0');
	my $public_experiments = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM experiment WHERE restricted=0 AND deleted=0');
	my $public_notebooks = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM list WHERE restricted=0 AND deleted=0');
	my $public_groups = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM user_group WHERE deleted=0');
	push @data, ['public', $public_genomes, $public_experiments, $public_notebooks, undef, $public_groups];
	
	my $owned_genomes = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM genome JOIN user_connector ON genome_id=child_id WHERE restricted=0 AND deleted=0 AND parent_type=5 AND child_type=2 AND role_id=2');
	my $owned_experiments = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM experiment JOIN user_connector ON experiment_id=child_id WHERE restricted=0 AND deleted=0 AND parent_type=5 AND child_type=3 AND role_id=2');
	my $owned_notebooks = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM list JOIN user_connector ON list_id=child_id WHERE restricted=0 AND deleted=0 AND parent_type=5 AND child_type=1 AND role_id=2');
	my $owned_groups = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM user_group JOIN user_connector ON user_group_id=child_id WHERE deleted=0 AND parent_type=5 AND child_type=6 AND user_connector.role_id=2');
	push @data, ['public (owned)', $owned_genomes, $owned_experiments, $owned_notebooks, undef, $owned_groups];
	
	my $shared_genomes = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM (SELECT 1 FROM genome JOIN user_connector ON genome_id=child_id WHERE deleted=0 AND restricted=1 AND parent_type=5 AND child_type=2 AND role_id!=2 GROUP BY genome_id) a');
	my $shared_experiments = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM (SELECT 1 FROM experiment JOIN user_connector ON experiment_id=child_id WHERE deleted=0 AND restricted=1 AND parent_type=5 AND child_type=3 AND role_id!=2 GROUP BY experiment_id) a');
	my $shared_notebooks = $db->storage->dbh->selectrow_array('SELECT FORMAT(COUNT(*), 0) FROM (SELECT 1 FROM list JOIN user_connector ON list_id=child_id WHERE deleted=0 AND restricted=1 AND parent_type=5 AND child_type=1 AND role_id!=2 GROUP BY list_id) a');
	push @data, ['restricted and shared', $shared_genomes, $shared_experiments, $shared_notebooks, undef, undef];

	return encode_json({
		data => \@data,
		bPaginate => 0,
	});
}

sub get_jobs_table {
	return encode_json({
		data => $db->storage->dbh->selectall_arrayref(
			"SELECT page,FORMAT(COUNT(*), 0) FROM log WHERE type != 0 AND page IN ('API','CoGeBlast','GEvo','LoadAnnotation','LoadExperiment','LoadGenome','SynFind','SynMap','SynMap2','SynMap3D') GROUP BY page"
		),
		bPaginate => 0
	});
}

sub add_row {
	my ($user_id, $user_name, $counts) = @_;
	my @jobs = ('API','CoGeBlast','GEvo','LoadAnnotation','LoadExperiment','LoadGenome','SynFind','SynMap','SynMap2','SynMap3D');
	my @row;
	push @row, '<a href="#" onclick="user_job_plot(' . $user_id . ',\'' . $user_name . '\')">' . $user_name . '</a>';
	foreach (@jobs) {
		push @row, $counts->{$_};
	}
	return \@row;
}

sub get_user_jobs_table {
	my @data;
	my %counts;
	my $public_jobs = $db->storage->dbh->selectall_arrayref(
		"SELECT page,FORMAT(COUNT(*), 0) FROM log WHERE type != 0 AND user_id = 0 AND page IN ('API','CoGeBlast','GEvo','LoadAnnotation','LoadExperiment','LoadGenome','SynFind','SynMap','SynMap2','SynMap3D') GROUP BY page"
	);
	foreach (@$public_jobs) {
		$counts{$_->[0]} = $_->[1];
	}
	push @data, add_row(0, ' public', \%counts);
	%counts = ();
	my $user_jobs = $db->storage->dbh->selectall_arrayref(
		"SELECT user.user_id,user_name,page,FORMAT(COUNT(*), 0) FROM user JOIN log ON user.user_id=log.user_id WHERE type != 0 AND page IN ('API','CoGeBlast','GEvo','LoadAnnotation','LoadExperiment','LoadGenome','SynFind','SynMap','SynMap2','SynMap3D') GROUP BY user_id,page"
	);
	my $user_id;
	my $user_name;
	foreach (@$user_jobs) {
		if ($user_id != $_->[0]) {
			push @data, add_row($user_id, $user_name, \%counts);
			$user_id = $_->[0];
			$user_name = $_->[1];
			%counts = ();
		}
		$counts{$_->[2]} = $_->[3];
	}
	push @data, add_row($user_id, $user_name, \%counts) if $user_name;

	return encode_json({
		data => \@data,
		bPaginate => 0
	});
}

sub get_user_jobs {
	my %opts = @_;
    my $user_id = $opts{user_id};
	my $user_jobs = $db->storage->dbh->selectall_arrayref("SELECT DATE(time),page,COUNT(*) FROM log WHERE user_id=" . $user_id . " AND type != 0 AND page IN ('API','CoGeBlast','GEvo','LoadAnnotation','LoadExperiment','LoadGenome','SynFind','SynMap','SynMap2','SynMap3D') GROUP BY DATE(time),page");
	my %jobs;
	foreach (@$user_jobs) {
		my $job = $jobs{$_->[1]};
		if ($job) {
			push @$job, [$_->[0], $_->[2]];
		}
		else {
			$jobs{$_->[1]} = [[$_->[0], $_->[2]]];
		}
	}
	my @data;
	foreach my $name (keys %jobs) {
		my $job = $jobs{$name};
		my @x;
		my @y;
		foreach my $run (@$job) {
			push @x, $run->[0];
			push @y, $run->[1];
		}
		push @data, { name => $name, x => \@x, y => \@y, type => 'bar' };
	}

	return encode_json(\@data);
}

####
#TAXONOMY STUFF

sub gen_tree_json {
	my %taxonomic_tree  = (
		name => "root", 
		children => [],
	);
	
	### 
	# Helper functions:
	###

	*check_tree = sub {
		my $search_term = shift;
		if (scalar(@_) == 0) {
			return undef;
		}
		
		my @new_roots;
		foreach my $root (@_) {
			foreach my $child (@{$root->{children}}) {
				push (@new_roots, $child);
				if(lc $search_term eq lc $child->{name}) {
					return $child;
				}
			}
		}
		return check_tree($search_term, @new_roots);
	};
	
	*gen_subtree = sub {
		my @array = @{$_[0]};
		if (scalar(@array) > 0) {
			my $hash = {
				name => $array[0],
				children => [],
			};
			shift @array;
			if (scalar(@array) > 0) {
				push ($hash->{children}, gen_subtree(\@array));
			}
			return $hash;
		}
	};

	*add_fix = sub {
		my $add_tree = $_[0];
		my $move_tree;
		my $i;
		for ($i = scalar(@{$taxonomic_tree{children}} - 1); $i > -1; $i--) {
			if ((lc $add_tree->{name}) eq (lc @{$taxonomic_tree{children}}[$i]->{name})) {
				$move_tree = splice(@{$taxonomic_tree{children}}, $i, 1);
				
				foreach my $child (@{$move_tree->{children}}) {
					add_to_tree($add_tree, $child);
				}
			}
		}
		#check for further inconsistencies
		foreach my $child (@{$add_tree->{children}}) {
			add_fix($child);
		}
	};
	
	*add_to_tree = sub {
		my $root = $_[0];
		my $sub_tree = $_[1];
	    my $root_children = \@{$root->{children}}; #array reference
		my $top_value = $sub_tree->{name};
		my $find = check_tree($top_value, $root);
		if (!$find) {
			#check for inconsistencies
			add_fix($sub_tree);
			
			#add the subtree as a child of the root
			push ($root_children, $sub_tree);
		} else {
			#recurse to the next level of the tree
			foreach my $child (@{$sub_tree->{children}}) {
				add_to_tree($find, $child);
			}
		}
	};
	
	###
	# Main code
	###

	my @organism_descriptions = CoGeDBI::get_table($db->storage->dbh, 'organism', undef, {description => "\"%;%\""}, {description => " like "});
	my $organism_descriptions = $organism_descriptions[0];
	
	foreach my $organism (values %$organism_descriptions) {
		my $taxon_string = $organism->{description};
		my @taxons = split(/\s*;\s*/, $taxon_string);
		
		if($taxons[0]) {
			$taxons[0] =~ s/^\s+|\s+$//g;
			my $size = scalar @taxons;
			$taxons[$size - 1] =~ s/^\s+|\s+$//g;
			my $sub_tree = gen_subtree(\@taxons);
			add_to_tree(\%taxonomic_tree, $sub_tree);
		}
	}
    
	return encode_json(\%taxonomic_tree);
}

	
####
#DATABASE TAB
	
sub get_total_queries {
	my $results = CoGeDBI::get_total_queries($db->storage->dbh);
	my $results2 = CoGeDBI::get_uptime($db->storage->dbh);	
    return encode_json({Queries => $results->{Queries}->{Value}, Uptime => $results2->{Uptime}->{Value}});
}
