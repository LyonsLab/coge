#! /usr/bin/perl -w

use strict;
use CGI;

#use CGI::Ajax;
use CoGeX;
use CoGe::Accessory::Web;
use HTML::Template;

no warnings 'redefine';

use vars qw($P $PAGE_TITLE $PAGE_NAME $USER $coge %FUNCTION $FORM $LINK);

$PAGE_TITLE = 'Groups';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

%FUNCTION = (
    get_groups_for_user      => \&get_groups_for_user,
    create_group             => \&create_group,
    delete_group             => \&delete_group,
    add_genome_to_group      => \&add_genome_to_group,
    remove_genome_from_group => \&remove_genome_from_group,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    #$template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER       => $name );
    $template->param( TITLE      => qq{Manage User Groups},
    				  PAGE_TITLE => $PAGE_TITLE,
    				  PAGE_LINK  => $LINK,
    				  HOME       => $P->{SERVER},
                      HELP       => 'Groups',
                      WIKI_URL   => $P->{WIKI_URL} || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( PAGE_NAME => $FORM->url );
    $template->param( MAIN      => 1 );
    my $groups = get_groups_for_user();
    $template->param( MAIN_TABLE => $groups );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

sub get_roles {
    my @roles;
    foreach my $role ( $coge->resultset('Role')->all() ) {
        next if $role->name =~ /admin/i && !$USER->is_admin;
        next if $role->name =~ /owner/i && !$USER->is_admin;
        my $name = $role->name;
        $name .= ": " . $role->description if $role->description;
        push @roles, { RID => $role->id, NAME => $name };
    }
    return [ reverse @roles ];    # put "reader" first
}

#sub add_genome_to_group {
#	my %opts  = @_;
#	my $ugid  = $opts{ugid};
#	my $dsgid = $opts{dsgid};
#	return "DSGID and/or UGID not specified" unless $dsgid && $ugid;
#	my $is_owner = $USER->is_owner( dsg => $dsgid );
#	$is_owner = 1 if $USER->is_admin;
#	return "You are not the owner of this genome" unless $is_owner;
#	my ($ugdc) =
#	  $coge->resultset('UserGroupDataConnector')
#	  ->create( { genome_id => $dsgid, user_group_id => $ugid } );
#	return 1;
#}
#
#sub remove_genome_from_group {
#	my %opts  = @_;
#	my $ugid  = $opts{ugid};
#	my $dsgid = $opts{dsgid};
#	return "DSGID and/or UGID not specified" unless $dsgid && $ugid;
#	my ($ugdc) =
#	  $coge->resultset('UserGroupDataConnector')
#	  ->search( { genome_id => $dsgid, user_group_id => $ugid } );
#	$ugdc->delete;
#	return 1;
#}

sub create_group {
    my %opts = @_;
    return "You need to be a registered user to create a user group!"
      unless $USER->id;
    my $name = $opts{name};
    return "No specified name!" unless $name;
    my $desc = $opts{desc};
    my $rid  = $opts{rid};

    unless ($rid) {
        my $role = $coge->resultset('Role')->find( { name => "Reader" } );
        $rid = $role->id;
    }

    # Create new group
    my $group = $coge->resultset('UserGroup')->create(
        {
            creator_user_id => $USER->id,
            name            => $name,
            description     => $desc,
            role_id         => $rid
        }
    );
    return unless $group;

# Make user the owner of new group
# my $conn = $coge->resultset('UserGroupConnector')->create( { user_id => $USER->id, user_group_id => $group->id } );
    my $conn = $coge->resultset('UserConnector')->create(
        {
            child_id    => $group->id,
            child_type  => 6,            #FIXME hardcoded to "group"
            parent_id   => $USER->id,
            parent_type => 5,            #FIXME hardcoded to "user"
            role_id     => 2             #FIXME hardcoded to "owner"
        }
    );
    return unless $conn;

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => 'create user group id' . $group->id
    );

    return 1;
}

sub delete_group {
    my %opts = @_;
    my $ugid = $opts{ugid};

    my $group = $coge->resultset('UserGroup')->find($ugid);
    return unless $group;
    if ( $group->locked && !$USER->is_admin ) {
        return "This is a locked group.  Admin permission is needed to delete.";
    }

    # OK, now delete the group
    $group->deleted(1);
    $group->update;

    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => "$PAGE_TITLE",
        description => 'delete user group id' . $group->id
    );

	return;
}

sub get_groups_for_user {
    my %opts = @_;

    my @group_list;
    if ( $USER->is_admin ) {
        @group_list = $coge->resultset('UserGroup')->all();
    }
    else {
        @group_list = $USER->groups;
    }

    my @groups;
    foreach my $group (@group_list) {

        #next if ($group->is_owner && !$USER->is_admin); # skip owner groups

        my $id          = $group->id;
        my $is_editable = $group->is_editable($USER);

        my %row;
        $row{NAME} =
            qq{<span class=link onclick='window.open("GroupView.pl?ugid=$id")'>}
          . $group->name
          . "</span>"
          . " (id$id)";
        $row{DESC} = $group->description if $group->description;
        my $role = $group->role->name;

        #		$role .= ": ".$group->role->description if $group->role->description;
        $row{ROLE} = $role;

        #		my $perm = join (", ", map {$_->name} $group->role->permissions);
        #		$row{PERM}=$perm;
        my @users;
        foreach
          my $user ( sort { $a->last_name cmp $b->last_name } $group->users )
        {
            my $display_name = $user->display_name;
            if ( $user->id == $group->creator_user_id )
            {    # is this user the creator?
                $display_name = '<b>'
                  . $display_name
                  . '</b> <span style="color: gray;">(creator)</span>';
            }
            push @users, $display_name;
        }

        #		push @users, "Self only" unless @users;
        $row{MEMBERS} = join( ",<br>", sort @users );
        my @data;
        foreach my $item (
            sort { $a->type->name cmp $b->type->name || $a->name cmp $b->name }
            $group->lists )
        {
            my $name =
              qq{<img src="picts/notebook-icon.png" width="15" height="15"/> }
              . $item->info_html;
            push @data, $name;
        }
        foreach my $item ( sort { $a->name cmp $b->name } $group->experiments )
        {
            my $name =
              qq{<img src="picts/testtube-icon.png" width="15" height="15"/> }
              . $item->info_html;
            push @data, $name;
        }
        foreach my $item ( sort { $a->name cmp $b->name } $group->genomes ) {
            my $name =
              qq{<img src="picts/dna-icon.png" width="15" height="15"/> }
              . $item->info_html;
            push @data, $name;
        }
        $row{DATA} = join( "<br>", @data );

        $row{BUTTONS} = 1;
        if ($is_editable) {
            $row{EDIT_BUTTON} =
"<span class='link ui-icon ui-icon-gear' onclick='window.open(\"GroupView.pl?ugid=$id\")'></span>";
            $row{DELETE_BUTTON} =
"<span class='link ui-icon ui-icon-trash' onclick=\"dialog_delete_group({ugid: "
              . $id
              . "});\"></span>";
        }
        else {
            $row{EDIT_BUTTON} = "<span class='ui-icon ui-icon-locked'></span>";
            $row{DELETE_BUTTON} =
              "<span class='ui-icon ui-icon-locked'></span>";
        }

        push @groups, \%row;
    }
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'Groups.tmpl' );
    $template->param( GROUP_TABLE => 1 );
    $template->param( GROUPS_LOOP => \@groups );
    $template->param( BUTTONS     => 1 );
    $template->param( ROLE_LOOP   => get_roles() );
    return $template->output;
}
