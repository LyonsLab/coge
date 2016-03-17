#! /usr/bin/perl -w

use strict;
use CGI;

use CoGeX;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON::XS;

no warnings 'redefine';

use vars
  qw( $P $PAGE_TITLE $PAGE_NAME $USER $coge %FUNCTION $FORM $MAX_SEARCH_RESULTS 
      $LINK $EMBED );

$PAGE_TITLE = 'GroupView';
$PAGE_NAME  = "$PAGE_TITLE.pl";
$MAX_SEARCH_RESULTS = 1000;

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    page_title => $PAGE_TITLE,
    cgi => $FORM
);

%FUNCTION = (
    get_group_info           => \&get_group_info,
    edit_group_info          => \&edit_group_info,
    update_group_info        => \&update_group_info,
    modify_users             => \&modify_users,
    add_user_to_group        => \&add_user_to_group,
    remove_user_from_group   => \&remove_user_from_group,
    delete_group             => \&delete_group,
    dialog_set_group_creator => \&dialog_set_group_creator,
    set_group_creator        => \&set_group_creator,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {  
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param( USER       => $USER->display_name || '',
                          TITLE      => qq{GroupView},
        			 	  PAGE_TITLE => qq{GroupView},
        				  PAGE_LINK  => $LINK,
        				  ADMIN_ONLY => $USER->is_admin,
        				  CAS_URL    => $P->{CAS_URL} || '',
        				  HOME       => $P->{SERVER},
                          HELP       => 'GroupView',
                          WIKI_URL   => $P->{WIKI_URL} || '' );
        $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    }
    
    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GroupView.tmpl' );
    $template->param( PAGE_NAME => $FORM->url );
    $template->param( MAIN      => 1 );
    my $ugid = $FORM->param('ugid');
    my %groups = map { $_->id => 1 } $USER->groups;

    unless ( $groups{$ugid} || $USER->is_admin ) {
        return "Group not found";
    }
    my ($group_info) = get_group_info( ugid => $ugid );
    $template->param( GROUP_INFO => $group_info );
    $template->param( UGID       => $ugid );

    $template->param( ADMIN_AREA    => $USER->is_admin );
    $template->param( ADMIN_BUTTONS => get_admin_functions() );

    #	my $open;
    #	$open = $FORM->param('open') if defined $FORM->param('open');
    #	my $box_open = $open ? 'true' : 'false';
    #	$template->param( EDIT_BOX_OPEN => $box_open );

    return $template->output;
}

sub get_admin_functions {
    my $ugid = $FORM->param('ugid');
    my $html =
qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="dialog_set_group_creator($ugid);">Set Group Creator</span>};
    return $html;
}

sub get_roles {
    my $current_role_id = shift;

    my @roles;
    foreach my $role ( $coge->resultset('Role')->all() ) {
        next if ( $role->name =~ /admin/i && !$USER->is_admin );
        next if ( $role->name =~ /owner/i && !$USER->is_admin );
        my $name = $role->name
          . ( $role->description ? ' (' . $role->description . ')' : '' );
        my $selected = '';
        $selected = 'selected="selected"' if ( $role->id == $current_role_id );
        push @roles,
          {
            ROLE_ID       => $role->id,
            ROLE_NAME     => $name,
            ROLE_SELECTED => $selected
          };
    }
    return \@roles;
}

sub get_group_info {
    my %opts = @_;
    my $ugid = $opts{ugid};
    return "Must have valid user group id.\n" unless ($ugid);

    my ($group) = $coge->resultset('UserGroup')->find($ugid);
    return "User group id$ugid does not exist.<br>" unless ($group);

    my $user_can_edit =
      ( $group->is_editable($USER)
          and ( not $group->locked or $USER->is_admin ) );

    my $html =
      $group->annotation_pretty_print_html;   #(allow_delete => $user_can_edit);
    if ($user_can_edit) {
        $html .= qq{<div class="coge-buttonset">};
        $html .=
qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="edit_group_info($ugid);">Edit Info</span>};
        $html .=
qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="modify_users($ugid);">Modify Users</span>};

#$html .= qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="add_lists({ugid: '$ugid'});">Add Notebook</span>};
        $html .=
			qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="delete_group();">} .
			($group->deleted ? 'Undelete' : 'Delete') . qq{</span>};
        $html .= qq{</div>};
    }

    return $html;
}

sub edit_group_info {
    my %opts = @_;
    my $ugid = $opts{ugid};
    return 0 unless $ugid;

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group && $group->is_editable($USER);

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GroupView.tmpl' );
    $template->param( EDIT_GROUP_INFO => 1 );
    $template->param( UGID            => $ugid );
    $template->param( NAME            => $group->name );
    $template->param( DESC            => $group->description );
    $template->param( ROLE_LOOP       => get_roles( $group->role->id ) );

    my %data;
    $data{title}  = 'Edit Group Info';
    $data{name}   = $group->name;
    $data{desc}   = $group->description;
    $data{output} = $template->output;

    return encode_json( \%data );
}

sub update_group_info {
    my %opts = @_;
    my $ugid = $opts{ugid};
    return 0 unless $ugid;
    my $name = $opts{name};
    return 0 unless $name;
    my $desc    = $opts{desc};
    my $role_id = $opts{roleid};
    return unless $role_id;

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return unless $group && $group->is_editable($USER);

    if ( $role_id != $group->role_id ) {
        foreach my $conn ( $group->child_connectors ) {
            $conn->role_id($role_id);
            $conn->update;
        }
        foreach my $conn ( $group->user_connectors ) {
            next if ( $conn->role->is_owner );
            $conn->role_id($role_id);
            $conn->update;
        }
    }

    $group->name($name);
    $group->description($desc) if $desc;
    $group->role_id($role_id);
    $group->update;

    return 1;
}

sub modify_users {
    my %opts = @_;
    my $ugid = $opts{ugid};
    return 0 unless $ugid;

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group && $group->is_editable($USER);

    my %data;
    $data{title} = 'Modify Users';
    my %users;
    my @users;

    foreach my $user (
        sort {
                 $a->last_name cmp $b->last_name
              || $a->user_name cmp $b->user_name
        } $group->users
      )
    {

        #		next if $user->user_name eq $USER->user_name; #skip self;
        push @users, { uid_name => $user->info, uid => $user->id };
        $users{ $user->id } = 1;
    }
    $data{users} = \@users;

    my @all_users;
    my $first = 1;
    foreach my $user (
        sort {
                 $a->last_name cmp $b->last_name
              || $a->user_name cmp $b->user_name
        } $coge->resultset('User')->all
      )
    {
        next if $users{ $user->id };    #skip users we already have
        push @all_users,
          {
            uid_name => $user->info,
            uid      => $user->id,
            selected => ( $first-- > 0 ? 'selected' : '' )
          };
    }
    $data{all_users} = \@all_users;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GroupView.tmpl' );
    $template->param( MODIFY_USERS  => 1 );
    $template->param( UGID          => $ugid );
    $template->param( UGID_LOOP     => $data{users} );
    $template->param( ALL_UGID_LOOP => $data{all_users} );

    $data{output} = $template->output;
    return encode_json( \%data );
}

sub add_user_to_group {
    my %opts = @_;
    my $ugid = $opts{ugid};
    my $uid  = $opts{uid};
    return "UGID and/or UID not specified" unless $ugid && $uid;

    #return 1 if $uid == $USER->id;

    # Find group and check permission to modify
    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group && $group->is_editable($USER);
    if ( $group->locked && !$USER->is_admin ) {
        return "This group is locked and cannot be modified.";
    }
    
    # Find user -- only needed for display name
    my $user = $coge->resultset('User')->find($uid);
    return unless $user;

    # Create user connection to group
    my $conn = $coge->resultset('UserConnector')->create(
        {
            parent_id   => $uid,
            parent_type => 5,                #FIXME hardcoded to "user"
            child_id    => $ugid,
            child_type  => 6,                #FIXME hardcoded to "group"
            role_id     => $group->role_id
        }
    );
    return 0 unless $conn;

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => 'added user ' . $user->info . ' to group ' . $group->info_html,
        parent_id   => $ugid,
        parent_type => 6 #FIXME magic number
    );

    return 1;
}

sub remove_user_from_group {
    my %opts = @_;
    my $ugid = $opts{ugid};
    my $uid  = $opts{uid};
    
    # Find user -- only needed for display name
    my $user = $coge->resultset('User')->find($uid);
    return unless $user;

    # Find group and check permission to modify
    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group && $group->is_editable($USER);
    if ( $group->locked && !$USER->is_admin ) {
        return "This group is locked and cannot be modified.";
    }
    if ( $uid == $USER->id && !$USER->is_admin ) {  # only allow this for admins
        return "Can't remove yourself from a group!";
    }
    if ( $group->owner && $uid == $group->owner->id && !$USER->is_admin )
    {                                               # only allow this for admins
        return "Can't remove the group owner!";
    }
    if ( $uid == $group->creator_user_id && !$USER->is_admin )
    {                                               # only allow this for admins
        return "Can't remove the group creator!";
    }
    
    # Remove all user connections to group
    foreach (
        $coge->resultset('UserConnector')->search(
            {
                parent_id   => $uid,
                parent_type => 5,                   #FIXME hardcoded to "user"
                child_id    => $ugid,
                child_type  => 6,                   #FIXME hardcoded to "group"
            }
        )
      )
    {
        $_->delete;
    }

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => 'removed user ' . $user->info . ' to group ' . $group->info_html,
        parent_id   => $ugid,
        parent_type => 6 #FIXME magic number
    );

    return 1;
}

sub delete_group {
    my %opts = @_;
    my $ugid = $opts{ugid};
    return "No UGID specified" unless $ugid;

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group && $group->is_editable($USER);

    if ( $group->locked && !$USER->is_admin ) {
        return "This is a locked group.  Admin permission is needed to delete.";
    }

    # OK, now delete the group
    $group->deleted(!$group->deleted); # do undelete if already deleted
    $group->update;

    return 1;
}

sub dialog_set_group_creator {
    my %opts = @_;
    my $ugid = $opts{ugid};
    return 0 unless $ugid;
    return 0 unless $USER->is_admin;

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group;

    my %data;
    $data{title} = 'Set Creator';
    my %users;
    my @users;

    my @all_users;
    foreach my $user (
        sort {
                 $a->last_name cmp $b->last_name
              || $a->user_name cmp $b->user_name
        } $coge->resultset('User')->all
      )
    {
        next if $users{ $user->id };    #skip users we already have
        my $name = $user->user_name;
        $name .= " : " . $user->first_name if $user->first_name;
        $name .= " " . $user->last_name
          if $user->first_name && $user->last_name;
        push @all_users, { uid_name => $name, uid => $user->id };
    }
    $data{all_users} = \@all_users;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GroupView.tmpl' );
    $template->param( SET_GROUP_CREATOR => 1 );
    $template->param( UGID              => $ugid );
    $template->param( ALL_UGID_LOOP     => $data{all_users} );

    $data{output} = $template->output;
    return encode_json( \%data );
}

sub set_group_creator {
    my %opts = @_;
    my $ugid = $opts{ugid};
    my $uid  = $opts{uid};
    return "UGID and/or UID not specified" unless $ugid && $uid;

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return 0 unless $group->is_editable($USER);

    if ( $group->locked && !$USER->is_admin ) {
        return "This is a locked group.  Admin permission is needed to modify.";
    }

    $group->creator_user_id($uid);
    $group->update();

    return 1;
}
