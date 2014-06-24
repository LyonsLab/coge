#! /usr/bin/perl -w

use strict;
use CGI;

use JSON::XS;
use HTML::Template;
use Sort::Versions;
use List::Util qw(first);
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;
use File::Path;
use File::stat;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Jex;
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;
use Benchmark;
no warnings 'redefine';

use vars qw($P $PAGE_TITLE $USER $LINK $coge %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS);

$PAGE_TITLE = 'User';

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
	page_title => $PAGE_TITLE,
    cgi => $FORM
);

# debug for fileupload:
# print STDERR $ENV{'REQUEST_METHOD'} . "\n" . $FORM->url . "\n" . Dumper($FORM->Vars) . "\n";	# debug
# print "data begin\n" . $FORM->param('POSTDATA') . "\ndata end\n" if ($FORM->param('POSTDATA'));

$MAX_SEARCH_RESULTS = 100;

my $node_types = CoGeX::node_types();

%ITEM_TYPE = (    # content/toc types
    all           		=> 100,
    mine          		=> 101,
    shared        		=> 102,
    activity      		=> 103,
    trash         		=> 104,
    activity_viz  		=> 105,
    activity_analyses	=> 106,
    user          		=> $node_types->{user},
    group         		=> $node_types->{group},
    notebook      		=> $node_types->{list},
    genome        		=> $node_types->{genome},
    experiment    		=> $node_types->{experiment}
);

%FUNCTION = (
    # get_logs				=> \&get_logs,
    upload_image_file               => \&upload_image_file,
    get_item_info                   => \&get_item_info,
    delete_items                    => \&delete_items,
    undelete_items                  => \&undelete_items,
    get_contents                    => \&get_contents,
    search_notebooks                => \&search_notebooks,
    add_items_to_notebook           => \&add_items_to_notebook,
    get_share_dialog                => \&get_share_dialog,
    search_share                    => \&search_share,
    add_items_to_user_or_group      => \&add_items_to_user_or_group,
    remove_items_from_user_or_group => \&remove_items_from_user_or_group,
	add_users_to_group				=> \&add_users_to_group,
	remove_user_from_group			=> \&remove_user_from_group,
    get_group_dialog                => \&get_group_dialog,
    change_group_role				=> \&change_group_role,
    send_items_to                   => \&send_items_to,
    create_new_group                => \&create_new_group,
    create_new_notebook             => \&create_new_notebook,
    toggle_star                     => \&toggle_star,
    cancel_jobs						=> \&cancel_jobs
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    #$template->param( TITLE      => 'User Profile' );
    $template->param( PAGE_TITLE => 'User Profile',
    				  PAGE_LINK  => $LINK,
    				  LOGO_PNG   => "$PAGE_TITLE-logo.png" );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    $template->param( ADJUST_BOX => 1 );

    return $template->output;
}

sub gen_body {
    if ( $USER->user_name eq 'public' ) {
        my $template =
          HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
        $template->param( LOGIN     => 1 );
        return $template->output;
    }

   # Other user specified as param, only allow access if collaborator
   # my $user = $USER;
   # my $uid = $FORM->param('uid');
   # if ($uid) {
   # 	return '' if (!$USER->is_admin && !$USER->has_collaborator($uid));
   # 	$user = $coge->resultset('User')->find($uid);
   # }
   # else {
   # 	my $uname = $FORM->param('name');
   # 	if ($uname) {
   # 		my $u = $coge->resultset('User')->find({user_name => $uname});
   # 		return '' if (!$u || (!$USER->is_admin && !$USER->has_collaborator($u)));
   # 		$user = $u;
   # 	}
   # }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        PAGE_NAME => "$PAGE_TITLE.pl",
        MAIN      => 1
    );

    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    $template->param(
        USER_NAME   => $USER->user_name,
        USER_ID     => $USER->id,
        FULL_NAME   => $USER->display_name,
        DESCRIPTION => $USER->description,
        EMAIL       => $USER->email,
        USER_IMAGE  => (
            $USER->image_id
            ? 'image.pl?id=' . $USER->image_id
            : 'picts/smiley_default.png'
        )
    );

    foreach ( keys %ITEM_TYPE ) {
        $template->param( 'ITEM_TYPE_' . uc($_) => $ITEM_TYPE{$_} );
    }

    # $template->param( LOGS => get_logs() );
    $template->param(
        TOC            => get_toc(),
        CONTENTS       => get_contents( html_only => 1 ),
        ROLES          => get_roles('reader'),
        NOTEBOOK_TYPES => get_notebook_types('mixed')
    );

    return $template->output;
}

sub get_item_info {
    my %opts      = @_;
    my $item_spec = $opts{item_spec};
    return unless $item_spec;
    my ( $item_id, $item_type ) = $item_spec =~ /content_(\d+)_(\d+)/;
    return unless ( $item_id and defined $item_type );
    my $timestamp = $opts{timestamp};

    # print STDERR "get_item_info: $item_id $item_type\n";

    my $html;
    if ( $item_type == $ITEM_TYPE{group} ) {
        my $group = $coge->resultset('UserGroup')->find($item_id);
        return unless $group;
        return unless ( $USER->is_admin or $group->has_member($USER) );

        $html .=
            '<b>Group id'
          . $group->id
          . '</b><br>'
          . '<b>Name:</b> '
          . $group->name . '<br>'
          . '<b>Description:</b> '
          . $group->description . '<br>'
          . '<b>Role:</b> '
          . $group->role->name . '<br>'
          . '<b>Members:</b><br>';
        foreach ( sort usercmp $group->users ) {
            $html .=
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
              . $_->display_name . ' ('
              . $_->user_name . ')' . '<br>';
        }
    }
    elsif ( $item_type == $ITEM_TYPE{notebook} ) {
        my $notebook = $coge->resultset('List')->find($item_id);
        return unless $USER->has_access_to_list($notebook);

        my $group_str = join( '<br>',
            sort map { $_->name } $USER->groups_with_access($notebook) );
        $html .=
            '<b>Notebook id'
          . $notebook->id
          . '</b><br>'
          . '<b>Name:</b> '
          . $notebook->name . '<br>'
          . '<b>Description:</b> '
          . $notebook->description . '<br>'
          . '<b>Contents:</b>'
          . '<div style="padding-left:20px;">'
          . $notebook->contents_summary_html
          . '</div>'
          . '<b>Groups with access:</b><br>'
          . '<div style="padding-left:20px;">'
          . ( $group_str ? $group_str : 'None' ) . '<br>'
          . '</div>'
          . '<b>Users with access:</b><br>'
          . '<div style="padding-left:20px;">';
        if ( $notebook->restricted ) {
            $html .= join( '<br>',
                sort map { $_->display_name . ' (' . $_->user_name . ')' }
                  $USER->users_with_access($notebook) );
        }
        else {
            $html .= 'Everyone';
        }
        $html .= '</div>';
    }
    elsif ( $item_type == $ITEM_TYPE{genome} ) {
        my $genome = $coge->resultset('Genome')->find($item_id);
        return unless ( $USER->has_access_to_genome($genome) );

        my $date =
          ( $genome->datasets ? $genome->datasets()->[0]->date : 'unknown' );
        my $group_str = join( '<br>',
            sort map { $_->name } $USER->groups_with_access($genome) );
        $html .=
            '<b>Genome id'
          . $genome->id
          . '</b><br>'
          . '<b>Organism: </b>'
          . $genome->organism->name . '<br>'
          . '<b>Name:</b> '
          . ( $genome->name ? $genome->name : '' ) . '<br>'
          . '<b>Description:</b> '
          . ( $genome->description ? $genome->description : '' ) . '<br>'
          . '<b>Version:</b> '
          . $genome->version . '<br>'
          . '<b>Type:</b> '
          . ( $genome->type ? $genome->type->name : '' ) . '<br>'
          . '<b>Source:</b> '
          . ( $genome->source ? $genome->source->[0]->name : '' ) . '<br>'
          . '<b>Created:</b> '
          . $date . '<br>'
          . '<b>Experiments:</b> '
          . $genome->experiments . '<br>'
          . '<b>Groups with access:</b><br>'
          . '<div style="padding-left:20px;">'
          . ( $group_str ? $group_str : 'None' ) . '<br>'
          . '</div>'
          . '<b>Users with access:</b><br>'
          . '<div style="padding-left:20px;">';
        if ( $genome->restricted ) {
            $html .= join( '<br>',
                sort map { $_->display_name . ' (' . $_->user_name . ')' }
                  $USER->users_with_access($genome) );
        }
        else {
            $html .= 'Everyone';
        }
        $html .= '</div>';
    }
    elsif ( $item_type == $ITEM_TYPE{experiment} ) {
        my $experiment = $coge->resultset('Experiment')->find($item_id);
        return unless $USER->has_access_to_experiment($experiment);

        my $group_str = join( '<br>',
            sort map { $_->name } $USER->groups_with_access($experiment) );
        $html .=
            '<b>Experiment id'
          . $experiment->id
          . '</b><br>'
          . '<b>Name:</b> '
          . $experiment->name . '<br>'
          . '<b>Description:</b> '
          . $experiment->description . '<br>'
          . '<b>Version:</b> '
          . $experiment->version . '<br>'
          . '<b>Source:</b> '
          . ( $experiment->source ? $experiment->source->name : '' ) . '<br>'
          . '<b>Created:</b> '
          . $experiment->date . '<br>'
          . '<b>Genome:</b> '
          . $experiment->genome->info . '<br>'
          . '<b>Groups with access:</b><br>'
          . '<div style="padding-left:20px;">'
          . ( $group_str ? $group_str : 'None' ) . '<br>'
          . '</div>'
          . '<b>Users with access:</b><br>'
          . '<div style="padding-left:20px;">';
        if ( $experiment->restricted ) {
            $html .= join( '<br>',
                sort map { $_->display_name . ' (' . $_->user_name . ')' }
                  $USER->users_with_access($experiment) );
        }
        else {
            $html .= 'Everyone';
        }
        $html .= '</div>';
    }

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub delete_items {
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
		my $type_name;

        #print STDERR "delete $item_id $item_type\n";
        if ( $item_type == $ITEM_TYPE{group} ) {
            my $group = $coge->resultset('UserGroup')->find($item_id);
		    return unless ( $group and $group->is_editable($USER) );
		    return if ( $group->locked and !$USER->is_admin );

			$group->deleted(1);
            $group->update;
			$type_name = 'user group';
        }
        elsif ( $item_type == $ITEM_TYPE{notebook} ) {
            my $notebook = $coge->resultset('List')->find($item_id);
            return unless $notebook;

            if ( !$notebook->locked
                and ( $USER->is_admin or $USER->is_owner( list => $notebook ) )
              )
            {
                $notebook->deleted(1);
                $notebook->update;
                $type_name = 'notebook';
            }
        }
        elsif ( $item_type == $ITEM_TYPE{genome} ) {
            my $genome = $coge->resultset('Genome')->find($item_id);
            return unless $genome;

            if ( $USER->is_admin or $USER->is_owner( dsg => $genome ) ) {
                $genome->deleted(1);
                $genome->update;
                $type_name = 'genome';
            }
        }
        elsif ( $item_type == $ITEM_TYPE{experiment} ) {
            my $experiment = $coge->resultset('Experiment')->find($item_id);
            return unless $experiment;

            if ( $USER->is_admin or $USER->is_owner( experiment => $experiment ) )
            {
                $experiment->deleted(1);
                $experiment->update;
                $type_name = 'experiment';
            }
        }
        
        # Record in log
        if ($type_name) {
			CoGe::Accessory::Web::log_history(
			    db          => $coge,
			    user_id     => $USER->id,
			    page        => $PAGE_TITLE,
			    description => "delete $type_name id$item_id"
			);
        }
    }
}

sub undelete_items {
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
        my $type_name;

        # print STDERR "undelete $item_id $item_type\n";
        if ( $item_type == $ITEM_TYPE{group} ) {
            my $group = $coge->resultset('UserGroup')->find($item_id);
            return unless $group;

            if ( $group->is_editable($USER) ) {
                $group->deleted(0);
                $group->update;
                $type_name = 'user group';
            }
        }
        elsif ( $item_type == $ITEM_TYPE{notebook} ) {
            my $notebook = $coge->resultset('List')->find($item_id);
            return unless $notebook;

            if ( $USER->is_admin or $USER->is_owner( list => $notebook ) ) {
                $notebook->deleted(0);
                $notebook->update;
                $type_name = 'notebook';
            }
        }                
        elsif ( $item_type == $ITEM_TYPE{genome} ) {
            my $genome = $coge->resultset('Genome')->find($item_id);
            return unless $genome;

            if ( $USER->is_admin or $USER->is_owner( dsg => $genome ) ) {
                $genome->deleted(0);
                $genome->update;
                $type_name = 'genome';
            }
        }
        elsif ( $item_type == $ITEM_TYPE{experiment} ) {
            my $experiment = $coge->resultset('Experiment')->find($item_id);
            return unless $experiment;

            if ( $USER->is_admin or $USER->is_owner( experiment => $experiment ) )
            {
                $experiment->deleted(0);
                $experiment->update;
                $type_name = 'experiment';
            }
        }
        
        # Record in log
        if ($type_name) {
			CoGe::Accessory::Web::log_history(
			    db          => $coge,
			    user_id     => $USER->id,
			    page        => $PAGE_TITLE,
			    description => "undelete $type_name id$item_id"
			);
        }
    }
}

sub cancel_jobs {
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

	my $jex = CoGe::Accessory::Jex->new( host => "localhost", port => 5151 );

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
        print STDERR "cancel $item_id $item_type\n";
        
        my $job = $coge->resultset('Job')->find($item_id);
	    if ( ( !$job || $job->user_id != $USER->id ) && !$USER->is_admin ) {
	        return;
	    }
	    
	    my $status = $jex->get_status( $job->id );
	    print STDERR "job " . $job->id . " status=$status\n";
	    if ( $status =~ /running/i ) {
	    	my $res = $jex->terminate( $job->id );
	    	if ( $res->{status} =~ /notfound/i ) {
	        	print STDERR "job " . $job->id . " termination error: status=" . $res->status . "\n";
	    	}
	    	else {
	    		$job->update( { status => 3 } );
	    	}
	        return encode_json( $res );
	    }
    }

	return 1;
}

sub get_roles {
    my $selected = shift;

    my $html;
    foreach my $role ( $coge->resultset('Role')->all() ) {
        next if $role->name =~ /admin/i && !$USER->is_admin;
        next if $role->name =~ /owner/i && !$USER->is_admin;
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

sub get_share_dialog {    #FIXME this routine needs to be optimized
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

    my ( %userconn, %notebooks );
    my $isPublic   = 0;
    my $isEditable = 1;
    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );

        # print STDERR "get_share $item_id $item_type\n";
        if ( $item_type == $ITEM_TYPE{genome} ) {
            my $genome = $coge->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );
            map { $userconn{ $_->parent_id } = $_ }
              ( $genome->user_connectors, $genome->group_connectors );
            map { $notebooks{ $_->id } = $_ } $genome->lists;
            $isPublic = 1 if ( not $genome->restricted );
            $isEditable = 0 if ( not $USER->is_owner_editor( dsg => $genome ) );
        }
        elsif ( $item_type == $ITEM_TYPE{experiment} ) {
            my $experiment = $coge->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);
            map { $userconn{ $_->id } = $_ }
              ( $experiment->user_connectors, $experiment->group_connectors );
            map { $notebooks{ $_->id } = $_ } $experiment->lists;
            $isPublic = 1 if ( not $experiment->restricted );
            $isEditable = 0
              if ( not $USER->is_owner_editor( experiment => $experiment ) );
        }
        elsif ( $item_type == $ITEM_TYPE{notebook} ) {
            my $notebook = $coge->resultset('List')->find($item_id);
            next unless $USER->has_access_to_list($notebook);
            map { $userconn{ $_->id } = $_ }
              ( $notebook->user_connectors, $notebook->group_connectors );
            $isPublic = 1 if ( not $notebook->restricted );
            $isEditable = 0
              if ( not $USER->is_owner_editor( list => $notebook ) );
        }
    }

    my ( %user_rows, %group_rows, %notebook_rows );
    foreach my $conn ( values %userconn ) {
        if ( $conn->is_parent_user ) {
            my $user = $conn->parent;
            $user_rows{ $user->id } = {
                USER_ITEM      => $user->id . ':' . $conn->parent_type,
                USER_FULL_NAME => $user->display_name,
                USER_NAME      => $user->name,
                USER_ROLE      => $conn->role->name,
                USER_DELETE    => $isEditable
                  && (!$conn->role->is_owner || $USER->is_admin)    # owner can't be removed unless admin
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
                GROUP_ITEM   => $group->id . ':' . $conn->parent_type,
                GROUP_NAME   => $group->name,
                GROUP_ROLE   => $group->role->name,
                GROUP_DELETE => $USER->is_owner_editor( group => $group->id ),
                GROUP_USER_LOOP => \@users
            };
        }
    }

    foreach my $notebook ( values %notebooks ) {
        my %users;

        foreach my $conn ( $notebook->user_connectors ) {
            my $user = $conn->parent;
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
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
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
        ROLES => get_roles('reader'),
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

    my ( %users, %roles, %creators, %owners, $lowest_role );
    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
        next unless ( $item_type == $ITEM_TYPE{group} ); # sanity check

        #print STDERR "get_group_dialog $item_id $item_type\n";
        my $group = $coge->resultset('UserGroup')->find($item_id);
        next unless ( $group and $group->is_editable($USER) );
		next if ( $group->locked and !$USER->is_admin );
		
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
            USER_ITEM      => $uid,
            USER_FULL_NAME => $user->display_name,
            USER_NAME      => $user->name,
            USER_ROLE      => ($role_name ? ' - ' . $role_name : ''),#$roles{$uid}->name),
            USER_DELETE    => !$owners{$uid} # owner can't be removed
        };
    }
    
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );

	# If no editable groups then show error dialog
    if (!$lowest_role) {
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
	    );
	}

    return $template->output;
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
          unless ( $_->user_name =~ /$search_term/i
            || $_->display_name =~ /$search_term/i );
        my $label = $_->display_name . ' (' . $_->user_name . ')';
        my $value = $_->id . ':' . $ITEM_TYPE{user};
        push @results, { 'label' => $label, 'value' => $value };
    }

    # Search for matching groups
    foreach ( $coge->resultset('UserGroup')->all ) {
        next unless ( $_->name =~ /$search_term/i );
        my $label = $_->name . ' (' . $_->role->name . ' group)';
        my $value = $_->id . ':' . $ITEM_TYPE{group};
        push @results, { 'label' => $label, 'value' => $value };
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

        # print STDERR "add_items_to_user_or_group $item_id $item_type\n";
        if ( $item_type == $ITEM_TYPE{genome} ) {
            my $genome = $coge->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );
            push @verified, $item;
        }
        elsif ( $item_type == $ITEM_TYPE{experiment} ) {
            my $experiment = $coge->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);
            push @verified, $item;
        }
        elsif ( $item_type == $ITEM_TYPE{notebook} ) {
            my $notebook = $coge->resultset('List')->find($item_id);
            next unless $USER->has_access_to_list($notebook);
            push @verified, $item;
        }
    }

    # Assign each item to user/group
    # print STDERR "add_items_to_user_or_group $target_id $target_type\n";
    #TODO verify that user can use specified role (for admin/owner roles)
    if ( $target_type == $ITEM_TYPE{user} ) {
        my $user = $coge->resultset('User')->find($target_id);
        return unless $user;

        foreach (@verified) {
            my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

            # print STDERR "   user: $item_id $item_type\n";

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
            return unless $conn;
        }
    }
    elsif ( $target_type == $ITEM_TYPE{group} ) {
        my $group = $coge->resultset('UserGroup')->find($target_id);
        return unless $group;

        foreach (@verified) {
            my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;

            # print STDERR "   group: $item_id $item_type\n";
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
    return unless $target_item;
    my $item_list = $opts{item_list};
    my @items = split( ',', $item_list );
    return unless @items;

    my ( $target_id, $target_type ) = $target_item =~ /(\d+)\:(\d+)/;
    next unless ( $target_id and $target_type );

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );

        # print STDERR "remove_item_from_user $item_id $item_type\n";
        if ( $item_type == $ITEM_TYPE{genome} ) {
            my $genome = $coge->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );

            my $conn = $coge->resultset('UserConnector')->find(
                {
                    parent_id   => $target_id,
                    parent_type => $target_type,
                    child_id    => $genome->id,
                    child_type  => $ITEM_TYPE{genome}
                }
            );
            return unless $conn;

            $conn->delete;
        }
        elsif ( $item_type == $ITEM_TYPE{experiment} ) {
            my $experiment = $coge->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);

            my $conn = $coge->resultset('UserConnector')->find(
                {
                    parent_id   => $target_id,
                    parent_type => $target_type,            #FIXME hardcoded
                    child_id    => $experiment->id,
                    child_type  => $ITEM_TYPE{experiment}
                }
            );
            return unless $conn;

            $conn->delete;
        }
        elsif ( $item_type == $ITEM_TYPE{notebook} ) {
            my $notebook = $coge->resultset('List')->find($item_id);
            next unless $USER->has_access_to_list($notebook);

            my $conn = $coge->resultset('UserConnector')->find(
                {
                    parent_id   => $target_id,
                    parent_type => $target_type,          #FIXME hardcoded
                    child_id    => $notebook->id,
                    child_type  => $ITEM_TYPE{notebook}
                }
            );
            return unless $conn;

            $conn->delete;
        }
    }

    return get_share_dialog( item_list => $item_list );
}

sub add_users_to_group {
    my %opts = @_;
    my @target_items = split( ',', $opts{target_items} );
    return unless @target_items;
    my $new_item = $opts{new_item};
    return unless $new_item;
    #print STDERR "add_users_to_group: $new_item @target_items\n";

	# Build a list of users to add to the target group
	my %users;
    my ( $item_id, $item_type ) = $new_item =~ /(\d+)\:(\d+)/;
    return unless ( $item_id and $item_type );

    if ( $item_type == $ITEM_TYPE{user} ) {
     	my $user = $coge->resultset('User')->find($item_id);
      	return unless $user;
       	$users{$user->id} = $user;
    }
    elsif ( $item_type == $ITEM_TYPE{group} ) {
      	my $group = $coge->resultset('UserGroup')->find($item_id);
        return unless $group;
        # TODO check that user has visibility of this group (one that they own or belong to)
		map { $users{$_->id} = $_ } $group->users;
     }
    
    # Add users to the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
		my ( $target_id, $target_type ) = $_ =~ /content_(\d+)_(\d+)/;
		#print STDERR "add_users_to_group $target_id\n";
	    next unless ( $target_id and $target_type );
	    next unless ( $target_type == $ITEM_TYPE{group} ); # sanity check
	    my $target_group = $coge->resultset('UserGroup')->find($target_id);
	    next unless ( $target_group and $target_group->is_editable($USER) );
	    next if ( $target_group->locked && !$USER->is_admin );

		# Add users to this target group
    	foreach my $user (values %users) {
    		# Check for existing user connection to target group
    		my $conn = $coge->resultset('UserConnector')->find(
		        {
		            parent_id   => $user->id,
		            parent_type => 5,                #FIXME hardcoded to "user"
		            child_id    => $target_id,
		            child_type  => 6                 #FIXME hardcoded to "group"
		        }
		    );
		    
		    # Create new user connection if one wasn't found
		    if (!$conn) {
			    $conn = $coge->resultset('UserConnector')->create(
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
		        db          => $coge,
		        user_id     => $USER->id,
		        page        => $PAGE_TITLE,
		        description => 'add user id' . $user->id . ' to group id' . $target_id
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
    #print STDERR "remove_user_from_group: $user_id @target_items\n";

	# Verify user
	my $user = $coge->resultset('User')->find($user_id);
    return unless $user;

    # Remove users from the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
		my ( $target_id, $target_type ) = $_ =~ /content_(\d+)_(\d+)/;
		#print STDERR "remove_user_from_group $target_id\n";
	    next unless ( $target_id and $target_type );
	    next unless ( $target_type == $ITEM_TYPE{group} ); # sanity check
	    my $target_group = $coge->resultset('UserGroup')->find($target_id);
	    next unless ( $target_group and $target_group->is_editable($USER) );
	    next if ( $target_group->locked && !$USER->is_admin );
	    
    	# Get user connection to target group
    	my $conn = $coge->resultset('UserConnector')->find(
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
		    db          => $coge,
		    user_id     => $USER->id,
		    page        => $PAGE_TITLE,
		    description => 'remove user id' . $user_id . ' from group id' . $target_id
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
    #print STDERR "change_group_role: $role_id @target_items\n";

	# Verify role
	my $role = $coge->resultset('Role')->find($role_id);
    return unless $role;
    
    # Change role for the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
		my ( $target_id, $target_type ) = $_ =~ /content_(\d+)_(\d+)/;
		#print STDERR "change_group_role $target_id\n";
	    next unless ( $target_id and $target_type );
	    next unless ( $target_type == $ITEM_TYPE{group} ); # sanity check
	    my $target_group = $coge->resultset('UserGroup')->find($target_id);
	    next unless ( $target_group and $target_group->is_editable($USER) );
	    next if ( $target_group->locked && !$USER->is_admin );

		$target_group->role_id($role_id);
		$target_group->update;
    }
    
	return get_group_dialog( item_list => $opts{target_items} );
}

sub send_items_to {
    my %opts      = @_;
    my $page_name = $opts{page_name};
    return unless $page_name;
    my $format    = $opts{format};
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

    my %fields;
    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
        push @{ $fields{$item_type} }, $item_id;
    }

    my $url;
    my $num = 1;
    foreach my $type ( keys %fields ) {
        my $name;
        if ( $type == $ITEM_TYPE{genome} ) {
            $name = 'dsgid';
        }
        elsif ( $type == $ITEM_TYPE{experiment} ) {
            $name = 'eid';
        }
        elsif ( $type == $ITEM_TYPE{notebook} ) {
            $name = 'nid';
        }

        if ( $format == 1 ) {    # numbered
            $url .= join( ';',
                map { $name . ( $num++ ) . '=' . $_ } @{ $fields{$type} } );
        }
        elsif ( $format == 2 ) {    # list
            $url .= $name . '=' . join( ',', @{ $fields{$type} } );
        }
        else {
            $url .= join( ';', map { $name . '=' . $_ } @{ $fields{$type} } );
        }
    }

    $url = $page_name . '.pl?' . $url if ($url);

    return $url;
}

sub get_toc {    # table of contents
    my @rows = (
    	{  	TOC_ITEM_ID       => $ITEM_TYPE{mine},
        	TOC_ITEM_INFO     => 'My Stuff',
        	TOC_ITEM_CHILDREN => 3
      	},
		{   TOC_ITEM_ID   => $ITEM_TYPE{notebook},
	        TOC_ITEM_INFO => 'Notebooks',
	        TOC_ITEM_ICON =>
	          '<img src="picts/notebook-icon.png" width="15" height="15"/>',
	        TOC_ITEM_INDENT => 20
      	},
      	{   TOC_ITEM_ID   => $ITEM_TYPE{genome},
	        TOC_ITEM_INFO => 'Genomes',
	        TOC_ITEM_ICON =>
	          '<img src="picts/dna-icon.png" width="15" height="15"/>',
	        TOC_ITEM_INDENT => 20
	    },
	    {   TOC_ITEM_ID   => $ITEM_TYPE{experiment},
	        TOC_ITEM_INFO => 'Experiments',
	        TOC_ITEM_ICON =>
	          '<img src="picts/testtube-icon.png" width="15" height="15"/>',
	        TOC_ITEM_INDENT => 20
	    },
    	{   TOC_ITEM_ID   => $ITEM_TYPE{shared},
	        TOC_ITEM_INFO => 'Shared with me'
	    },
		{ 	TOC_ITEM_ID => $ITEM_TYPE{group},
 	 		TOC_ITEM_INFO => 'Groups',
 			#TOC_ITEM_ICON => '<img src="picts/group-icon.png" width="15" height="15"/>' 
		},
    	{   TOC_ITEM_ID       => $ITEM_TYPE{activity},
	        TOC_ITEM_INFO     => 'Activity',
	        TOC_ITEM_CHILDREN => 2
	    },
	    {   TOC_ITEM_ID     => $ITEM_TYPE{activity_analyses},
	        TOC_ITEM_INFO   => 'Analyses',
	        TOC_ITEM_INDENT => 20
	    },
    	{   TOC_ITEM_ID     => $ITEM_TYPE{activity_viz},
	        TOC_ITEM_INFO   => 'Graph',
	        TOC_ITEM_INDENT => 20
	    },
    	{	TOC_ITEM_ID   => $ITEM_TYPE{trash},
        	TOC_ITEM_INFO => 'Trash'
      	}
    );

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( DO_TOC        => 1 );
    $template->param( TOC_ITEM_LOOP => \@rows );
    return $template->output;
}

sub get_contents {
    my %opts = @_;
    my $type = $opts{type};
    $type = $ITEM_TYPE{all} unless $type;
    my $last_update = $opts{last_update};
    $last_update = 0 if ( not defined $last_update );
    my $timestamp = $opts{timestamp};
    my $html_only = $opts{html_only};
    #print STDERR "get_contents $type $html_only $last_update\n";

    #use Time::HiRes qw ( time );
    #my $start_time = time;
    my @rows;

    # Get current time (according to database)
    my $update_time = $coge->storage->dbh_do(
        sub {
            my ( $storage, $dbh, @args ) = @_;
            $dbh->selectrow_array('SELECT NOW()+0');
        }
    );

    # Preload stuff speed (needed for genome/experiment sorting and info routines)
    #FIXME which would be faster, children_by_type_role_id or joins?
    my %sourceIdToName =
      map { $_->id => $_->name } $coge->resultset('DataSource')->all();

    # Get all items for this user (genomes, experiments, notebooks)
    #my $t1    = new Benchmark;
    my ( $children, $roles ) = $USER->children_by_type_role_id;
    #my $t2    = new Benchmark;
    #my $time = timestr( timediff( $t2, $t1 ) );
    #print STDERR "User.pl: Number of children\n";
    #foreach my $k (sort keys %$children) {
	#	print STDERR $k,"\t", scalar values $children->{$k},"\n";
    #}
    #$children={};
    #print STDERR "Time to run user->children_by_type_role_id: $time\n";
    #print STDERR "get_contents: time1=" . ((time - $start_time)*1000) . "\n";

	if (   $type == $ITEM_TYPE{all} 
		or $type == $ITEM_TYPE{group} ) 
	{
	 	foreach my $group ( sort {$a->name cmp $b->name} values %{ $children->{6} } ) { #FIXME hardcoded type
	 		push @rows, { CONTENTS_ITEM_ID => $group->id,
	 					  CONTENTS_ITEM_TYPE => $ITEM_TYPE{group},
	 					  CONTENTS_ITEM_DELETED => $group->deleted,
	 					  CONTENTS_ITEM_INFO => $group->info,
	 				  	  CONTENTS_ITEM_ICON => '<img src="picts/group-icon.png" width="15" height="15" style="vertical-align:middle;"/>',
	 				  	  CONTENTS_ITEM_LINK =>  'GroupView.pl?ugid=' . $group->id,
	 				  	  CONTENTS_ITEM_SELECTABLE => 1
	 		};
	 	}
	}

    if (   $type == $ITEM_TYPE{notebook}
        or $type == $ITEM_TYPE{all}
        or $type == $ITEM_TYPE{mine} )
    {
        foreach my $list ( sort listcmp values %{ $children->{1} } )
        {    #FIXME hardcoded type
            push @rows, {
                CONTENTS_ITEM_ID   => $list->id,
                CONTENTS_ITEM_TYPE => $ITEM_TYPE{notebook},
                CONTENTS_ITEM_DELETED => $list->deleted,
                CONTENTS_ITEM_SHARED =>
                  !$roles->{2}{ $list->id },    #FIXME hardcoded role id
                CONTENTS_ITEM_INFO => $list->info,
                CONTENTS_ITEM_ICON =>
'<img src="picts/notebook-icon.png" width="15" height="15" style="vertical-align:middle;"/>',
                CONTENTS_ITEM_LINK       => 'NotebookView.pl?nid=' . $list->id,
                CONTENTS_ITEM_SELECTABLE => 1
            };
        }
    }

    #print STDERR "get_contents: time2=" . ((time - $start_time)*1000) . "\n";
    if (   $type == $ITEM_TYPE{genome}
        or $type == $ITEM_TYPE{all}
        or $type == $ITEM_TYPE{mine} )
    {
        foreach my $genome ( sort genomecmp values %{ $children->{2} } )
        {    #FIXME hardcoded type
            push @rows, {
                CONTENTS_ITEM_ID      => $genome->id,
                CONTENTS_ITEM_TYPE    => $ITEM_TYPE{genome},
                CONTENTS_ITEM_DELETED => $genome->deleted,
                CONTENTS_ITEM_SHARED =>
                  !$roles->{2}{ $genome->id },    #FIXME hardcoded role id
                CONTENTS_ITEM_INFO => $genome->info,
                CONTENTS_ITEM_ICON =>
'<img src="picts/dna-icon.png" width="15" height="15" style="vertical-align:middle;"/>',
                CONTENTS_ITEM_LINK       => 'GenomeInfo.pl?gid=' . $genome->id,
                CONTENTS_ITEM_SELECTABLE => 1
            };
        }
    }

    #print STDERR "get_contents: time3=" . ((time - $start_time)*1000) . "\n";
    if (   $type == $ITEM_TYPE{experiment}
        or $type == $ITEM_TYPE{all}
        or $type == $ITEM_TYPE{mine} )
    {
        foreach my $experiment ( sort experimentcmp values %{ $children->{3} } )
        {    #FIXME hardcoded type
            push @rows, {
                CONTENTS_ITEM_ID      => $experiment->id,
                CONTENTS_ITEM_TYPE    => $ITEM_TYPE{experiment},
                CONTENTS_ITEM_DELETED => $experiment->deleted,
                CONTENTS_ITEM_SHARED =>
                  !$roles->{2}{ $experiment->id },    #FIXME hardcoded role id
                CONTENTS_ITEM_INFO => $experiment->info(
                    source => $sourceIdToName{ $experiment->data_source_id }
                ),
                CONTENTS_ITEM_ICON =>
'<img src="picts/testtube-icon.png" width="15" height="15" style="vertical-align:middle;"/>',
                CONTENTS_ITEM_LINK => 'ExperimentView.pl?eid='
                  . $experiment->id,
                CONTENTS_ITEM_SELECTABLE => 1
            };
        }
    }

    #print STDERR "get_contents: time4=" . ((time - $start_time)*1000) . "\n";
    if ( $type == $ITEM_TYPE{all} or $type == $ITEM_TYPE{activity} ) {
        foreach my $entry (
            $USER->logs(
                {
                    time => { '>=' => $last_update },
                    type => { '!=' => 0 }
                },    #FIXME hardcoded type
                { order_by => { -desc => 'time' } }
            )
          )
        {
            my $icon = '<img id="' . $entry->id . '" ' . ( $entry->is_important ? 'src="picts/star-full.png"' : 'src="picts/star-hollow.png"' ) . ' ' . 'width="15" height="15" style="vertical-align:middle;" ' . 'onclick="toggle_star(this);"' . '/>';
            push @rows,
              {
                CONTENTS_ITEM_ID   => $entry->id,
                CONTENTS_ITEM_TYPE => $ITEM_TYPE{activity},
                CONTENTS_ITEM_INFO => $entry->short_info,
                CONTENTS_ITEM_ICON => $icon,
                CONTENTS_ITEM_LINK => $entry->link
              };
        }
    }
    
    #print STDERR "get_contents: time5=" . ((time - $start_time)*1000) . "\n";
    if ( $type == $ITEM_TYPE{all} or $type == $ITEM_TYPE{activity_analyses} ) {
        push @rows, 
          { CONTENTS_ITEM_ID   => 0,
            CONTENTS_ITEM_TYPE => $ITEM_TYPE{activity_analyses},
            CONTENTS_ITEM_INFO => '<span class="small alert">Please note: This area is undergoing changes and the analyses shown do not reflect recent activity.</span>',
            #CONTENTS_ITEM_LINK => undef,
            CONTENTS_ITEM_SELECTABLE => 0
          };
        
        foreach my $entry (
            $USER->jobs(
                {},
                {
                	time => { '>=' => $last_update },
                	join => 'log',
                	prefetch => 'log',
                	order_by => { -desc => 'start_time' } 
                }
            )
          )
        {
            push @rows,
              {
                CONTENTS_ITEM_ID   => $entry->id,
                CONTENTS_ITEM_TYPE => $ITEM_TYPE{activity_analyses},
                CONTENTS_ITEM_INFO => $entry->info_html,
                CONTENTS_ITEM_LINK => $entry->link,
                CONTENTS_ITEM_SELECTABLE => 1
              };
        }
    }
    #print STDERR "get_contents: time6=" . ((time - $start_time)*1000) . "\n";
    
    if ($html_only) {    # only do this for initial page load, not polling
        my $user_id = $USER->id;
        my $job_list = 'cogeblast/synmap/gevo/synfind/loadgenome/loadexperiment/organismview/user';
        push @rows,
          {
            CONTENTS_ITEM_ID   => 0,
            CONTENTS_ITEM_TYPE => $ITEM_TYPE{activity_viz},
            CONTENTS_ITEM_INFO => qq{<iframe frameborder="0" width="100%" height="100%" scrolling="no" src="//genomevolution.org/blacktea/standalone/$user_id/$job_list"></iframe>} # FIXME: hardcoded server name
          };
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( DO_CONTENTS        => 1 );
    $template->param( CONTENTS_ITEM_LOOP => \@rows );
    my $html = $template->output;

    #print STDERR "get_contents: time=" . ((time - $start_time)*1000) . "\n";

    return $html if ($html_only);
    return encode_json(
        { timestamp => $timestamp, lastUpdate => $update_time, html => $html }
    );
}

# sub get_logs {
# 	return if ($USER->user_name eq "public");

# 	my %opts = @_;
# 	my $type = $opts{type};

# 	my @logs;
# 	if (!$type or $type eq 'recent') {
# 		@logs = $coge->resultset('Log')->search( { user_id => $USER->id, description => { 'not like' => 'page access' } }, { order_by => { -desc => 'time' } } ); # $user->logs;
# 		#my @logs = reverse $coge->resultset('Log')->search_literal( 'user_id = ' . $user->id . ' AND time >= DATE_SUB(NOW(), INTERVAL 1 HOUR)' );
# 	}
# 	else {
# 		@logs = $coge->resultset('Log')->search( { user_id => $USER->id, status => 1 }, { order_by => { -desc => 'time' } } );
# 	}

# 	my @rows;
# 	foreach (splice(@logs, 0, 100)) {
# 		push @rows, { LOG_TIME => $_->time,
# 					  LOG_PAGE => $_->page,
# 					  LOG_DESC => $_->description,
# 					  LOG_LINK => $_->link
# 					};
# 	}
# 	return if (not @rows);

# 	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
# 	$template->param( LOG_TABLE => 1 );
# 	$template->param( LOG_LOOP => \@rows );
# 	return $template->output;
# }

sub upload_image_file {
    return if ( $USER->user_name eq "public" );

    my %opts           = @_;
    my $image_filename = '' . $FORM->param('input_upload_file');
    my $fh             = $FORM->upload('input_upload_file');
    return if ( -s $fh > 2 * 1024 * 1024 );    # limit to 2MB

    #TODO delete old image

    # Create the image
    my $image;
    if ($fh) {

        #print STDERR "$image_filename size=" . (-s $fh) . "\n";
        read( $fh, my $contents, -s $fh );
        $image = $coge->resultset('Image')->create(
            {
                filename => $image_filename,
                image    => $contents
            }
        );
        return unless $image;

        # Link to user
        $USER->image_id( $image->id );
        $USER->update;
        return encode_json( { link => 'image.pl?id=' . $image->id } );
    }

    return;
}

sub search_notebooks
{    # FIXME this coded is dup'ed in CoGeBlast.pl and NotebookView.pl
    my %opts = @_;
    return if ( $USER->user_name eq 'public' );
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$search_term $timestamp\n";

    my @notebooks;
    my $num_results;
    my $group_str = join( ',', map { $_->id } $USER->groups );

    # Try to get all items if blank search term
    if ( !$search_term ) {
        my $sql = "locked=0"
          ;    # AND restricted=0 OR user_group_id IN ( $group_str ))"; # FIXME
        $num_results = $coge->resultset("List")->count_literal($sql);
        if ( $num_results < $MAX_SEARCH_RESULTS ) {
            foreach
              my $notebook ( $coge->resultset("List")->search_literal($sql) )
            {
                next unless $USER->has_access_to_list($notebook);
                push @notebooks, $notebook;
            }
        }
    }

    # Perform search
    else {

        # Get public lists and user's private lists
        $search_term = '%' . $search_term . '%';
        foreach my $notebook (
            $coge->resultset("List")->search_literal(
"locked=0 AND (name LIKE '$search_term' OR description LIKE '$search_term')"
            )
          )
        {
            next unless $USER->has_access_to_list($notebook);
            push @notebooks, $notebook;
        }
        $num_results = @notebooks;
    }

    # Limit number of results display
    if ( $num_results > $MAX_SEARCH_RESULTS ) {
        return encode_json(
            {
                timestamp => $timestamp,
                html =>
"<option>$num_results matches, please refine your search.</option>"
            }
        );
    }

    # Build select items out of results
    my $html;
    foreach my $n ( sort listcmp @notebooks ) {
        $html .=
          "<option value='" . $n->id . "'>" . $n->info . "</option><br>\n";
    }
    $html = "<option disabled='disabled'>No matches</option>" unless $html;

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub add_items_to_notebook {
    my %opts = @_;
    my $nid  = $opts{nid};
    return unless $nid;
    my $item_list = $opts{item_list};
    my @items = split( ',', $item_list );
    return unless @items;

    # print STDERR "add_items_to_notebook $nid $item_list\n";

    my $notebook = $coge->resultset('List')->find($nid);
    return unless $USER->has_access_to_list($notebook);

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /content_(\d+)_(\d+)/;
        next unless ( $item_id and $item_type );
        next
          unless ( $item_type eq $ITEM_TYPE{notebook}
            or $item_type eq $ITEM_TYPE{genome}
            or $item_type eq $ITEM_TYPE{experiment} );

        #TODO check access permission on each item

        # print STDERR "add_item_to_notebook $item_id $item_type\n";

        my $conn = $coge->resultset('ListConnector')->find_or_create(
            {
                parent_id  => $nid,
                child_id   => $item_id,
                child_type => $item_type
            }
        );
        return unless $conn;
    }

    return 1;
}

sub create_new_group {
    my %opts    = @_;
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    my $role_id = $opts{role_id};
    return unless $name && $role_id;

    return if ( $USER->user_name eq "public" );

    my $role = $coge->resultset('Role')->find($role_id);
    return unless $role;

    # Create the new group
    my $group = $coge->resultset('UserGroup')->create(
        {
            creator_user_id => $USER->id,
            name            => $name,
            description     => $desc,
            role_id         => $role->id
        }
    );
    return unless $group;

    # Set user as owner
    my $conn = $coge->resultset('UserConnector')->create(
        {
            parent_id   => $USER->id,
            parent_type => 5,            #FIXME hardcoded to "user"
            child_id    => $group->id,
            child_type  => 6,            #FIXME hardcoded to "group"
            role_id     => 2             #FIXME hardcoded to "owner"
        }
    );
    return unless $conn;

    # Add selected items to new group
    #TODO

    # Record in log
    CoGe::Accessory::Web::log_history(
		db          => $coge,
		user_id     => $USER->id,
		page        => $PAGE_TITLE,
		description => 'create user group id' . $group->id
	);

    return 1;
}

sub create_new_notebook {
    my %opts    = @_;
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    my $type_id = $opts{type_id};
    return unless $name && $type_id;
    my $item_list = $opts{item_list};    # optional
    return if ( $USER->user_name eq "public" );

    # Create the new list
    my $list = $coge->resultset('List')->create(
        {
            name         => $name,
            description  => $desc,
            list_type_id => $type_id,

            # user_group_id => $owner->id,
            restricted => 1
        }
    );
    return unless $list;

    # Set user as owner
    my $conn = $coge->resultset('UserConnector')->create(
        {
            parent_id   => $USER->id,
            parent_type => 5,           #FIXME hardcoded to "user"
            child_id    => $list->id,
            child_type  => 1,           #FIXME hardcoded to "list"
            role_id     => 2            #FIXME hardcoded to "owner"
        }
    );
    return unless $conn;

    # Add selected items to new notebook
    add_items_to_notebook( nid => $list->id, item_list => $item_list )
      if ($item_list);

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => "$PAGE_TITLE",
        description => 'create notebook id' . $list->id
    );

    return 1;
}

sub get_notebook_types {
    my $selected = shift;
    my $html;
    foreach my $type ( $coge->resultset('ListType')->all() ) {
        next
          if ( $type->name =~ /owner/i )
          ;    # reserve this type for system-created lists
        my $name =
          $type
          ->name;    # . ($type->description ? ": " . $type->description : '');
        $html .=
            '<option value="'
          . $type->id . '" '
          . (    $type->id eq $selected
              || $type->name =~ /$selected/i ? 'selected' : '' )
          . '>'
          . $name
          . '</option>';
    }
    return $html;
}

sub toggle_star {
    my %opts   = @_;
    my $log_id = $opts{log_id};

    my $entry = $coge->resultset('Log')->find($log_id);
    return '' unless $entry;

    my $status = $entry->status;
    $entry->status( not $status );
    $entry->update();

    return not $status;
}

# FIXME these comparison routines are duplicated elsewhere
sub genomecmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->organism->name cmp $b->organism->name
      || versioncmp( $b->version, $a->version )
      || $a->type->id <=> $b->type->id
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

sub experimentcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    versioncmp( $b->version, $a->version )
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

sub listcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->name cmp $b->name;
}

sub groupcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->name cmp $b->name;
}

sub usercmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->display_name cmp $b->display_name;
}
