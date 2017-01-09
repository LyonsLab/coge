#! /usr/bin/perl -w

use v5.10;
use strict;

use Switch;
use CGI;
use Data::Validate::URI qw(is_uri);
use JSON::XS;
use HTML::Template;
use Sort::Versions;
use Time::Piece;
use List::Util qw(first);
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;
use Benchmark;
use File::Path;
use File::stat;
use File::Spec::Functions qw(catdir catfile);
use CoGeDBI;
use CoGeX;
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;
use CoGe::Accessory::Utils qw(format_time_diff js_escape html_escape commify);
use CoGe::Accessory::Web;
use CoGe::JEX::Jex;
use CoGe::Core::Item qw(get_item);
use CoGe::Core::Notebook qw(notebookcmp);
use CoGe::Core::Experiment qw(experimentcmp);
use CoGe::Core::Genome qw(genomecmp);
use CoGe::Core::Metadata qw(create_annotations create_image);
use CoGe::Core::Favorites;
no warnings 'redefine';

use vars qw(
    $CONF $PAGE_TITLE $USER $LINK $DB %FUNCTION $FORM %ITEM_TYPE
    $MAX_SEARCH_RESULTS $JEX $node_types
);

$PAGE_TITLE = 'User';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    page_title => $PAGE_TITLE,
    cgi => $FORM
);

# Admins have ability to simulate other users using the "user_id" query parameter
my $user_id = $FORM->Vars->{'user_id'};
if (defined $user_id && $USER->is_admin && $user_id != $USER->id) {
    my $user = $DB->resultset('User')->find($user_id);
    if (defined $user) {
        print STDERR "Switching to user '", $user->name, "'\n";
        $USER = $user;
    }
}

$JEX = CoGe::JEX::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );

$MAX_SEARCH_RESULTS = 100;

$node_types = CoGeX::node_types();

# Content types
# note: adding new values is okay, but don't change existing values or 
# will break running User pages in the field
%ITEM_TYPE = (
    all           		=> 100,
    mine          		=> 101,
    shared        		=> 102,
    activity_summary	=> 103,
    trash         		=> 104,
    activity_viz  		=> 105,
    activity_analyses	=> 106,
    activity_loads      => 107,
    metadata			=> 108,
    user          		=> $node_types->{user},
    group         		=> $node_types->{group},
    notebook      		=> $node_types->{list},
    genome        		=> $node_types->{genome},
    experiment    		=> $node_types->{experiment}
);

%FUNCTION = (
    upload_image_file               => \&upload_image_file,
    get_item_info                   => \&get_item_info,
    favorite_items                  => \&favorite_items,
    delete_items                    => \&delete_items,
    undelete_items                  => \&undelete_items,
    get_contents                    => \&get_contents,
    search_notebooks                => \&search_notebooks,
    add_items_to_notebook           => \&add_items_to_notebook,
    get_share_dialog                => \&get_share_dialog,
    search_share                    => \&search_share,
    make_items_public               => \&make_items_public,
    add_items_to_user_or_group      => \&add_items_to_user_or_group,
    remove_items_from_user_or_group => \&remove_items_from_user_or_group,
    add_users_to_group			    => \&add_users_to_group,
    remove_user_from_group		    => \&remove_user_from_group,
    get_group_dialog                => \&get_group_dialog,
    change_group_role			    => \&change_group_role,
    send_items_to                   => \&send_items_to,
    create_new_group                => \&create_new_group,
    create_new_notebook             => \&create_new_notebook,
    toggle_star                     => \&toggle_star,
    cancel_job				        => \&cancel_job,
    comment_job                     => \&comment_job,
    upload_metadata					=> \&upload_metadata
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( 
        USER       => $USER->display_name || '',
        PAGE_TITLE => 'My Data',
	    TITLE      => "My Data",
        PAGE_LINK  => $LINK,
		SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
		HOME       => $CONF->{SERVER},
        HELP       => 'User',
        WIKI_URL   => $CONF->{WIKI_URL} || '',
        ADMIN_ONLY => $USER->is_admin,
        CAS_URL    => $CONF->{CAS_URL} || '',
        COOKIE_NAME => $CONF->{COOKIE_NAME} || ''
    );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    return $template->output;
}

sub gen_body {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    
    if ( $USER->is_public ) {
        $template->param( PAGE_NAME => "$PAGE_TITLE.pl",
                          LOGIN     => 1 );
        return $template->output;
    }
    
    foreach ( keys %ITEM_TYPE ) {
        $template->param( 'ITEM_TYPE.' . uc($_) => $ITEM_TYPE{$_} );
    }

    $template->param(
        PAGE_NAME   => "$PAGE_TITLE.pl",
        MAIN        => 1,
        USER_NAME   => $USER->user_name,
        USER_ID     => $USER->id,
        FULL_NAME   => $USER->display_name || '',
        EMAIL       => $USER->email,
        ITEM_TYPE     => encode_json(\%ITEM_TYPE),
        ROLES          => get_roles('reader'),
        #NOTEBOOK_TYPES => get_notebook_types('mixed') # mdb removed 12/14/16 COGE-800
    );

    return $template->output;
}

sub get_item_info { #TODO move into API and render on client-side -- this is also called by search-results.js
    my %opts      = @_;
    my $item_id = $opts{item_id};
    my $item_type = $opts{item_type};
    return unless ($item_id and $item_type);
    my $timestamp = $opts{timestamp};
    #print STDERR "get_item_info: ", Dumper \%opts, "\n";

    my $html;
    if ( $item_type eq 'group' ) {
        my $group = $DB->resultset('UserGroup')->find($item_id);
        return unless $group;
        return unless ( $USER->is_admin or $group->has_member($USER) );
        
        my $creation = ($group->creator_user_id ? $group->creator->display_name  . ' ' : '') . ($group->date ne '0000-00-00 00:00:00' ? $group->date : '');

        $html =
            '<b>Group id' . $group->id . '</b><br>'
          . '<b>Name:</b> ' . $group->name . '<br>'
          . '<b>Description:</b> ' . $group->description . '<br>'
          . '<b>Role:</b> ' . $group->role->name . '<br>'
          . ($creation ? '<b>Creator:</b> ' . $creation . '<br>' : '')
          . '<b>Members:</b><br>';
        foreach ( sort usercmp $group->users ) {
            $html .=
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
              . $_->display_name . ' (' . $_->user_name . ')' . '<br>';
        }

        my $view_link = qq{open_item('GroupView.pl?ugid=$item_id');};
        
        $html .= qq{<div><b>Tools:</b><br>}
            . qq{<div style="padding-left:20px;">}
            . qq{<span class="link" onclick="$view_link" title="View contents">View group</span><br>};
        $html .= qq{<span class="link" onclick="group_dialog();" title="Edit membership">Edit membership</span><br>}
            if ( ( $group->is_editable($USER) and ( not $group->locked or $USER->is_admin ) ) );
        $html .= qq{</div></div>};
    }
    elsif ( $item_type eq 'notebook' ) {
        my $notebook = $DB->resultset('List')->find($item_id);
        return unless $USER->has_access_to_list($notebook);

        my $group_str = join( '<br>',
            sort map { $_->name } $USER->groups_with_access($notebook) );
        $html = '<b>Notebook id' . $notebook->id . '</b><br>';
        $html .= '<b>Name:</b> ' . $notebook->name . '<br>';
        $html .= '<b>Description:</b> ' . $notebook->description . '<br>' if $notebook->description;
        $html .= '<b>Contents:</b><div style="padding-left:20px;">';
        $html .= $notebook->contents_summary_html if $notebook->contents_summary_html;
        $html .= '</div><b>Groups with access:</b><br><div style="padding-left:20px;">';
        $html .= ( $group_str ? $group_str : 'None' ) . '<br>';
        $html .= '</div><b>Users with access:</b><br><div style="padding-left:20px;">';
        if ( $notebook->restricted ) {
            $html .= join( '<br>',
                sort map { $_->display_name . ' (' . $_->user_name . ')' }
                  $USER->users_with_access($notebook) );
        }
        else {
            $html .= 'Everyone';
        }
        $html .= '</div>';
        
        my $view_link = qq{open_item('NotebookView.pl?lid=$item_id');};

        $html .= qq{<div><b>Tools:</b><br>}
            . qq{<div style="padding-left:20px;">}
            . qq{<span class="link" onclick="$view_link" title="View contents">View details</span><br>};

        $html .= qq{<span class="link" onclick="share_dialog();" title="Share with other users or user groups">Share</span><br>}
            if ($notebook->is_editable($USER));

        $html .= qq{</div></div>};
    }
    elsif ( $item_type eq 'organism' ) {
        my $organism = $DB->resultset('Organism')->find($item_id);
        return unless $organism;

        $html =
            '<b>Organism id' . $organism->id . '</b><br>'
          . '<b>Name:</b> ' . $organism->name . '<br>'
          . '<b>Description:</b> ' . $organism->description . '<br>'
          . '<b>Public genomes:</b> ' . ($organism->public_genomes > 0 ? $organism->public_genomes : 'none') . '<br>';

        my $view_link = qq{open_item('OrganismView.pl?oid=$item_id');};

        $html .= qq{<div><b>Tools:</b><br>}
            . qq{<div style="padding-left:20px;">}
            . qq{<span class="link" onclick="$view_link" title="View/edit metadata and contents">View details</span><br>}
            . qq{</div></div>};
    }
    elsif ( $item_type eq 'genome' || $item_type eq 'favorite' ) {
        my $genome = $DB->resultset('Genome')->find($item_id);
        return unless ( $USER->has_access_to_genome($genome) );

        my $date = ( $genome->datasets ? $genome->datasets()->[0]->date : 'unknown' );
        my $group_str = join( '<br>', sort map { $_->name } $USER->groups_with_access($genome) );
        $html =
            '<b>Genome id' . $genome->id . '</b><br>'
          . '<b>Organism: </b>' . $genome->organism->name . '<br>'
          . '<b>Chromosomes: </b>' . commify($genome->chromosome_count) . '<br>'
          . '<b>Name:</b> ' . ( $genome->name ? $genome->name : '' ) . '<br>'
          . '<b>Description:</b> ' . ( $genome->description ? $genome->description : '' ) . '<br>'
          . '<b>Version:</b> ' . $genome->version . '<br>'
          . '<b>Type:</b> ' . ( $genome->type ? $genome->type->name : '' ) . '<br>'
          . '<b>Source:</b> ' . ( $genome->source ? $genome->source->[0]->name : '' ) . '<br>'
          . '<b>Created:</b> ' . $date . '<br>'
          . '<b>Annotated:</b> ' . ($genome->has_gene_features ? 'yes' : 'no') . '<br>'
          . '<b>Experiments:</b> ' . $genome->experiments . '<br>'
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
        
        my $edit_link = qq{open_item('GenomeInfo.pl?gid=$item_id');};
        my $view_link = qq{open_item('GenomeView.pl?gid=$item_id&tracks=sequence%2Cfeatures');};
        my $load_link = qq{open_item('LoadAnnotation.pl?gid=$item_id');};
        
        $html .= qq{<div><b>Tools:</b><br>}
            . qq{<div style="padding-left:20px;">}
            . qq{<span class="link" onclick="$edit_link" title="View/edit metadata, annotation, and experiments">View details</span><br>}
            . qq{<span class="link" onclick="$view_link" title="Browse sequence, annotation, and experiment tracks">Browse</span><br>};
        if ($genome->is_editable($USER)) {
            $html .= qq{<span class="link" onclick="share_dialog();" title="Share with other users or user groups">Share</span><br>}
                . qq{<span class="link" onclick="$load_link" title="Load gene annotation">Load annotation</span><br>};
        }
        $html .= qq{<span class="link" onclick="add_to_notebook_dialog();" title="Add to a notebook">Add to notebook</span><br>}
            if (!$USER->is_public);
        $html .= qq{</div></div>};
    }
    elsif ( $item_type eq 'feature' ) {
        my $feature = $DB->resultset('Feature')->find($item_id);
        return unless $feature;

        $html =
            '<b>Feature id' . $feature->id . '</b><br>'
          . '<b>Type:</b> ' . join(', ', map { $_->name } $feature->type) . '<br>'
          . '<b>Primary name:</b> ' . ($feature->primary_name ? join (', ', map { $_->name } $feature->primary_name) : 'none') . '<br>'
          . '<b>Other names:</b> ' . join(', ', $feature->names) . '<br>'
          . '<b>Genome:</b> ' . join(', ', map { $_->info } $feature->genomes) . '<br>'
          . '<b>Location:</b> ' . $feature->chromosome . ':' . $feature->start . '-' . $feature->stop . '<br>'
          . '<b>Length:</b> ' . abs($feature->stop - $feature->start + 1) . ' bp';

        my $view_link = qq{open_item('FeatView.pl?fid=$item_id');};

        $html .= qq{<div><b>Tools:</b><br>}
            . qq{<div style="padding-left:20px;">}
            . qq{<span class="link" onclick="$view_link" title="View details">View details</span><br>}
            . qq{</div></div>};
    }
    elsif ( $item_type eq 'experiment' ) {
        my $experiment = $DB->resultset('Experiment')->find($item_id);
        return unless $USER->has_access_to_experiment($experiment);

        my $group_str = join( '<br>', sort map { $_->name } $USER->groups_with_access($experiment) );
        $html =
            '<b>Experiment id' . $experiment->id . '</b><br>'
          . '<b>Name:</b> ' . $experiment->name . '<br>'
          . '<b>Description:</b> ' . $experiment->description . '<br>'
          . '<b>Version:</b> ' . $experiment->version . '<br>'
          . '<b>Source:</b> ' . ( $experiment->source ? $experiment->source->name : '' ) . '<br>'
          . '<b>Created:</b> ' . $experiment->date . '<br>'
          . '<b>Genome:</b> ' . $experiment->genome->info . '<br>'
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
        
        my $gid = $experiment->genome_id;
        my $edit_link = qq{open_item('ExperimentView.pl?eid=$item_id');};
        my $view_link = qq{open_item('GenomeView.pl?gid=$gid&tracks=experiment$item_id');};
        
        $html .= qq{<div><b>Tools:</b><br>}
            . qq{<div style="padding-left:20px;">}
            . qq{<span class="link" onclick="$edit_link" title="View/edit metadata">View details</span><br>}
            . qq{<span class="link" onclick="$view_link" title="Browse track data">Browse data</span><br>};
        if ($experiment->is_editable($USER)) {
            $html .= qq{<span class="link" onclick="share_dialog();" title="Share with other users or user groups">Share</span><br>};
        }
        $html .= qq{<span class="link" onclick="add_to_notebook_dialog();" title="Add a notebook">Add to notebook</span><br>}
            if (!$USER->is_public);
        $html .= qq{</div></div>};
    }
    elsif ( $item_type eq 'analyses' or $item_type eq 'loads' ) {
        my $log = $DB->resultset('Log')->find($item_id);
        return unless $log;
        
        $html =
            '<b>Workflow id' . $log->parent_id . '</b><br>'
          . '<b>Type:</b> ' . $log->page . '<br>'
          . '<b>Description:</b> ' . $log->description . '<br>'
          . '<b>Date:</b> ' . $log->time . '<br>';
        $html .= '<b>Comment:</b> ' . $log->comment . '<br>' if $log->comment;
        $html .= qq{<div><b>Tools:</b><br>}
          . qq{<div style="padding-left:20px;">}
          . qq{<a href="} . $log->link . qq{" target=_blank>Open result</a>}
          . qq{</div>}
    }
    else {
        $html = '';
    }

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub favorite_items {
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );

        my $item_obj = get_item($DB, $item_id, $item_type);
        next unless $item_obj;
        
        my $favorites = CoGe::Core::Favorites->new(user => $USER);
        my $is_favorited = $favorites->toggle($item_obj);

        # Record in log
        my $item_type_code = $node_types->{$item_type};
        CoGe::Accessory::Web::log_history(
            db          => $DB,
            user_id     => $USER->id,
            page        => $PAGE_TITLE,
            description => ($is_favorited ? 'Favorited' : 'Unfavorited') . " $item_type " . $item_obj->info_html,
            parent_id   => $item_id,
            parent_type => $item_type_code
        );
    }
}

sub delete_items {
    my %opts      = @_;
    my $item_list = $opts{item_list};
    my @items     = split( ',', $item_list );
    return unless @items;

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );
		my $type_name;

        #print STDERR "delete $item_id $item_type\n";
        my $item_obj;
        if ( $item_type eq 'group' ) {
            $item_obj = $DB->resultset('UserGroup')->find($item_id);
		    return unless ( $item_obj and $item_obj->is_editable($USER) );
		    return if ( $item_obj->locked and !$USER->is_admin );

			$item_obj->deleted(1);
            $item_obj->update;
			$type_name = 'user group';
        }
        elsif ( $item_type eq 'notebook' ) {
            $item_obj = $DB->resultset('List')->find($item_id);
            return unless $item_obj;

            if ( !$item_obj->locked
                and ( $USER->is_admin or $USER->is_owner( list => $item_obj ) )
              )
            {
                $item_obj->deleted(1);
                $item_obj->update;
                $type_name = 'notebook';
            }
        }
        elsif ( $item_type eq 'genome' ) {
            $item_obj = $DB->resultset('Genome')->find($item_id);
            return unless $item_obj;

            if ( $USER->is_admin or $USER->is_owner( dsg => $item_obj ) ) {
                $item_obj->deleted(1);
                $item_obj->update;
                $type_name = 'genome';
            }
        }
        elsif ( $item_type eq 'experiment' ) {
            $item_obj = $DB->resultset('Experiment')->find($item_id);
            return unless $item_obj;

            if ( $USER->is_admin or $USER->is_owner( experiment => $item_obj ) )
            {
                $item_obj->deleted(1);
                $item_obj->update;
                $type_name = 'experiment';
            }
        }

        # Record in log
        if ($type_name) {
            my $item_type_code = $node_types->{$item_type};
			CoGe::Accessory::Web::log_history(
			    db          => $DB,
			    user_id     => $USER->id,
			    page        => $PAGE_TITLE,
			    description => "deleted $type_name " . $item_obj->info_html,
			    parent_id   => $item_id,
			    parent_type => $item_type_code
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
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );
        my $type_name;

        #print STDERR "undelete $item_id $item_type\n";
        my $item_obj;
        if ( $item_type eq 'group' ) {
            $item_obj = $DB->resultset('UserGroup')->find($item_id);
            return unless $item_obj;

            if ( $item_obj->is_editable($USER) ) {
                $item_obj->deleted(0);
                $item_obj->update;
                $type_name = 'user group';
            }
        }
        elsif ( $item_type eq 'notebook' ) {
            $item_obj = $DB->resultset('List')->find($item_id);
            return unless $item_obj;

            if ( $USER->is_admin or $USER->is_owner( list => $item_obj ) ) {
                $item_obj->deleted(0);
                $item_obj->update;
                $type_name = 'notebook';
            }
        }
        elsif ( $item_type eq 'genome' ) {
            $item_obj = $DB->resultset('Genome')->find($item_id);
            return unless $item_obj;

            if ( $USER->is_admin or $USER->is_owner( dsg => $item_obj ) ) {
                $item_obj->deleted(0);
                $item_obj->update;
                $type_name = 'genome';
            }
        }
        elsif ( $item_type eq 'experiment' ) {
            $item_obj = $DB->resultset('Experiment')->find($item_id);
            return unless $item_obj;

            if ( $USER->is_admin or $USER->is_owner( experiment => $item_obj ) )
            {
                $item_obj->deleted(0);
                $item_obj->update;
                $type_name = 'experiment';
            }
        }
        else {
            print STDERR "User::undelete_items error: unknown type\n";
        }

        # Record in log
        if ($type_name) {
            my $item_type_code = $node_types->{$item_type};
			CoGe::Accessory::Web::log_history(
			    db          => $DB,
			    user_id     => $USER->id,
			    page        => $PAGE_TITLE,
			    description => "undeleted $type_name " . $item_obj->info_html,
			    parent_id   => $item_id,
			    parent_type => $item_type_code
			);
        }
    }
}

sub cancel_job {
    my %opts      = @_;
    my $workflow_id = $opts{workflow_id};
    return 0 unless $workflow_id;

    my $status = $JEX->get_status( $workflow_id );
    if ( $status =~ /running/i ) {
        my $res = $JEX->terminate( $workflow_id );
        if ( $res =~ /notfound/i ) {
            print STDERR "job " . $workflow_id . " termination error: status=" . $res . "\n";
        }
        return 0;
    }

    return 1;
}

sub comment_job {
    my %opts      = @_;
    my $log_id = $opts{log_id};
    my $comment = $opts{comment};
    return 0 unless $log_id;
    
    warn "comment_job: $log_id $comment";

    my $log = $USER->logs->find($log_id);
    return 0 unless $log;

    $log->comment($comment);
    $log->update;

    return 1;
}

sub get_roles {
    my $selected = shift;

    my $html;
    foreach my $role ( $DB->resultset('Role')->all() ) {
        next if ($role->name =~ /admin/i);# && !$USER->is_admin); # skip admin
        next if ($role->name =~ /owner/i && !$USER->is_admin); # skip owner
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
    #print STDERR 'get_share_dialog: ', $item_list, "\n";
    my @items     = split( ',', $item_list );
    return unless @items;

    my ( %userconn, %notebooks );
    my $isPublic   = 0;
    my $isEditable = 1;
    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );

        # print STDERR "get_share $item_id $item_type\n";
        if ( $item_type eq 'genome' ) {
            my $genome = $DB->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );
            map { $userconn{ $_->parent_id } = $_ }
              ( $genome->user_connectors, $genome->group_connectors );
            map { $notebooks{ $_->id } = $_ } $genome->lists;
            $isPublic = 1 if ( not $genome->restricted );
            $isEditable = 0 if ( not $USER->is_owner_editor( dsg => $genome ) );
        }
        elsif ( $item_type eq 'experiment' ) {
            my $experiment = $DB->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);
            map { $userconn{ $_->id } = $_ }
              ( $experiment->user_connectors, $experiment->group_connectors );
            map { $notebooks{ $_->id } = $_ } $experiment->lists;
            $isPublic = 1 if ( not $experiment->restricted );
            $isEditable = 0
              if ( not $USER->is_owner_editor( experiment => $experiment ) );
        }
        elsif ( $item_type eq 'notebook' ) {
            my $notebook = $DB->resultset('List')->find($item_id);
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
	   next unless $conn->role; #EL added 10/21/14 so solve a problem for a user where the share dialog wouldn't appear for his genomes, and an fatal error was being thrown due to role->name with role being undefined.  Genomes with problems were IDs: 24576, 24721, 24518, 24515, 24564, 24566, 24568, 24562, 24571 
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

            if (not $group->deleted) {
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
    }

    foreach my $notebook ( values %notebooks ) {
        next if $notebook->deleted;
        
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

    my $template = HTML::Template->new( filename => catfile($CONF->{TMPLDIR}, 'partials', 'user_common.tmpl') );
    $template->param(
        SHARE_DIALOG => 1,
        IS_EDITABLE  => $USER->is_admin || $isEditable,
        IS_RESTRICTED => !$isPublic,
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
        $template->param( ACCESS_MSG => 'Everyone (publicly available)' );
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
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );
        next unless ( $item_type eq 'group' ); # sanity check

        #print STDERR "get_group_dialog $item_id $item_type\n";
        my $group = $DB->resultset('UserGroup')->find($item_id);
        next unless ( $group and $group->is_editable($USER) );
		next if ( $group->locked and !$USER->is_admin );

		my $role = $group->role;
		$lowest_role = $role if (!$lowest_role or $role->is_lower($lowest_role));
		my $creator_id = $group->creator_user_id;
		my $owner_id;
        $owner_id = $group->owner->id if $group->owner;
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

    my $template = HTML::Template->new( filename => catfile($CONF->{TMPLDIR}, 'partials', 'user_common.tmpl') );

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
    foreach ( $DB->resultset('User')->all ) {
        next
          unless ($search_term && ( escape($_->user_name) =~ /$search_term/i
            || escape($_->display_name) =~ /$search_term/i ));
        next if ($_->user_name eq $USER->user_name); # mdb added 5/11/16 -- prevent user from sharing with themselves
        my $label = $_->display_name . ' (' . $_->user_name . ')';
        my $value = $_->id . ':' . 'user';
        push @results, { 'label' => $label, 'value' => $value };
    }

    # Search for matching groups
    foreach ( $DB->resultset('UserGroup')->all ) {
        next unless ( $search_term && $_->name && escape($_->name) =~ /$search_term/i );
        my $label = $_->name . ' (' . $_->role->name . ' group)';
        my $value = $_->id . ':' . 'group';
        push @results, { 'label' => $label, 'value' => $value };
    }

    return encode_json( { timestamp => $timestamp, items => \@results } );
}

sub make_items_public {
    my %opts        = @_;
    my $make_public = $opts{make_public};
    $make_public = 1 unless defined $make_public;
    my $item_list = $opts{item_list};
    my @items = split( ',', $item_list );
    return unless @items;
    #print STDERR "make_items_public ", $item_list, " ", $make_public, "\n";
    
    # Verify that user has access to each item then make public
    foreach my $item (@items) {
        my ( $item_id, $item_type ) = $item =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );

        #print STDERR "make_items_public $item_id $item_type\n";
        if ( $item_type eq 'genome' ) {
            my $genome = $DB->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );
            $genome->restricted(!$make_public);
            $genome->update();
        }
        elsif ( $item_type eq 'experiment' ) {
            my $experiment = $DB->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);
            $experiment->restricted(!$make_public);
            $experiment->update();
        }
        elsif ( $item_type eq 'notebook' ) {
            my $notebook = $DB->resultset('List')->find($item_id);
            next unless $USER->has_access_to_list($notebook);
            $notebook->restricted(!$make_public);
            $notebook->update();
        }
    }
    
    return get_share_dialog( item_list => $item_list );
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
    #print STDERR "add_items_to_user_or_group ", $target_item, " ", $item_list, "\n";

    my ( $target_id, $target_type ) = $target_item =~ /(\d+)\:(\w+)/;
    return unless ( $target_id and $target_type );

    # Verify that user has access to each item
    my @verified;
    foreach my $item (@items) {
        my ( $item_id, $item_type ) = $item =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );

        # print STDERR "add_items_to_user_or_group $item_id $item_type\n";
        if ( $item_type eq 'genome' ) {
            my $genome = $DB->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );
            push @verified, { id => $item_id, type => $ITEM_TYPE{genome}, type_name => $item_type, info => $genome->info_html };
        }
        elsif ( $item_type eq 'experiment' ) {
            my $experiment = $DB->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);
            push @verified, { id => $item_id, type => $ITEM_TYPE{experiment}, type_name => $item_type, info => $experiment->info_html };
        }
        elsif ( $item_type eq 'notebook' ) {
            my $notebook = $DB->resultset('List')->find($item_id);
            next unless $USER->has_access_to_list($notebook);
            push @verified, { id => $item_id, type => $ITEM_TYPE{notebook}, type_name => $item_type, info => $notebook->info_html };
        }
    }

    # Assign each item to user/group
    # print STDERR "add_items_to_user_or_group $target_id $target_type\n";
    #TODO verify that user can use specified role (for admin/owner roles)
    if ( $target_type eq 'user' ) { #$ITEM_TYPE{user} ) {
        my $user = $DB->resultset('User')->find($target_id);
        return unless $user;
        next if ($user->user_name eq $USER->user_name); # mdb added 5/11/16 -- prevent user from sharing with themselves

        foreach (@verified) {
            my ( $item_id, $item_type, $item_type_name, $item_info ) = ( $_->{id}, $_->{type}, $_->{type_name}, $_->{info} );
            # print STDERR "   user: $item_id $item_type\n";

            # Remove previous connection
            foreach ($DB->resultset('UserConnector')->search(
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
            my $conn = $DB->resultset('UserConnector')->create(
                {
                    parent_id   => $target_id,
                    parent_type => 5,            # FIXME hardcoded
                    child_id    => $item_id,
                    child_type  => $item_type,
                    role_id     => $role_id
                }
            );
            return unless $conn;
            
            # Record in log
            CoGe::Accessory::Web::log_history(
                db          => $DB,
                user_id     => $USER->id,
                page        => $PAGE_TITLE,
                description => "shared (role=$role_id) $item_type_name $item_info with user " . $user->info,
                parent_id   => $target_id,
                parent_type => 6 #FIXME magic number
            );
        }
    }
    elsif ( $target_type eq 'group' ) {
        my $group = $DB->resultset('UserGroup')->find($target_id);
        return unless $group;

        foreach (@verified) {
            my ( $item_id, $item_type, $item_type_name, $item_info ) = ( $_->{id}, $_->{type}, $_->{type_name}, $_->{info} );

            # print STDERR "   group: $item_id $item_type\n";
            my $conn = $DB->resultset('UserConnector')->find_or_create(
                {
                    parent_id   => $target_id,
                    parent_type => 6,                # FIXME hardcoded
                    child_id    => $item_id,
                    child_type  => $item_type,
                    role_id     => $group->role_id
                }
            );
            return unless $conn;
            
            # Record in log
            CoGe::Accessory::Web::log_history(
                db          => $DB,
                user_id     => $USER->id,
                page        => $PAGE_TITLE,
                description => "shared $item_type_name $item_info with group " . $group->info_html,
                parent_id   => $target_id,
                parent_type => 6 #FIXME magic number
            );
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
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );

        # print STDERR "remove_item_from_user $item_id $item_type\n";
        if ( $item_type eq 'genome' ) {
            my $genome = $DB->resultset('Genome')->find($item_id);
            next unless ( $USER->has_access_to_genome($genome) );

            my $conn = $DB->resultset('UserConnector')->find(
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
        elsif ( $item_type eq 'experiment' ) {
            my $experiment = $DB->resultset('Experiment')->find($item_id);
            next unless $USER->has_access_to_experiment($experiment);

            my $conn = $DB->resultset('UserConnector')->find(
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
        elsif ( $item_type eq 'notebook' ) {
            my $notebook = $DB->resultset('List')->find($item_id);
            next unless $USER->has_access_to_list($notebook);

            my $conn = $DB->resultset('UserConnector')->find(
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
    my ( $item_id, $item_type ) = $new_item =~ /(\d+)\:(\w+)/;
    return unless ( $item_id and $item_type );

    if ( $item_type eq 'user' ) {
     	my $user = $DB->resultset('User')->find($item_id);
      	return unless $user;
       	$users{$user->id} = $user;
    }
    elsif ( $item_type eq 'group' ) {
      	my $group = $DB->resultset('UserGroup')->find($item_id);
        return unless $group;
        # TODO check that user has visibility of this group (one that they own or belong to)
		map { $users{$_->id} = $_ } $group->users;
     }

    # Add users to the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
		my ( $target_id, $target_type ) = $_ =~ /(\d+)_(\w+)/;
		#print STDERR "add_users_to_group $target_id\n";
	    next unless ( $target_id and $target_type );
	    next unless ( $target_type eq 'group' ); # sanity check
	    my $target_group = $DB->resultset('UserGroup')->find($target_id);
	    next unless ( $target_group and $target_group->is_editable($USER) );
	    next if ( $target_group->locked && !$USER->is_admin );

		# Add users to this target group
    	foreach my $user (values %users) {
    		# Check for existing user connection to target group
    		my $conn = $DB->resultset('UserConnector')->find(
		        {
		            parent_id   => $user->id,
		            parent_type => 5,                #FIXME hardcoded to "user"
		            child_id    => $target_id,
		            child_type  => 6                 #FIXME hardcoded to "group"
		        }
		    );

		    # Create new user connection if one wasn't found
		    if (!$conn) {
			    $conn = $DB->resultset('UserConnector')->create(
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
		        db          => $DB,
		        user_id     => $USER->id,
		        page        => $PAGE_TITLE,
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
    #print STDERR "remove_user_from_group: $user_id @target_items\n";

	# Verify user
	my $user = $DB->resultset('User')->find($user_id);
    return unless $user;

    # Remove users from the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
		my ( $target_id, $target_type ) = $_ =~ /(\d+)_(\w+)/;
		#print STDERR "remove_user_from_group $target_id\n";
	    next unless ( $target_id and $target_type );
	    next unless ( $target_type eq 'group' ); # sanity check
	    my $target_group = $DB->resultset('UserGroup')->find($target_id);
	    next unless ( $target_group and $target_group->is_editable($USER) );
	    next if ( $target_group->locked && !$USER->is_admin );

    	# Get user connection to target group
    	my $conn = $DB->resultset('UserConnector')->find(
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
		    db          => $DB,
		    user_id     => $USER->id,
		    page        => $PAGE_TITLE,
		    description => 'removed user ' . $user->info . ' from group ' . $target_group->info_html,
		    parent_id   => $target_id,
		    parent_type => 6 #FIXME magic number
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
	my $role = $DB->resultset('Role')->find($role_id);
    return unless $role;

    # Change role for the target groups
    foreach (@target_items) {
        # Find target group and check permission to modify
		my ( $target_id, $target_type ) = $_ =~ /(\d+)_(\w+)/;
		#print STDERR "change_group_role $target_id\n";
	    next unless ( $target_id and $target_type );
	    next unless ( $target_type eq 'group' ); # sanity check
	    my $target_group = $DB->resultset('UserGroup')->find($target_id);
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
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );
        push @{ $fields{$item_type} }, $item_id;
    }

    my $url;
    my $num = 1;
    foreach my $type ( keys %fields ) {
        my $name;
        if ( $type eq 'genome' ) {
            $name = 'gid';
        }
        elsif ( $type eq 'experiment' ) {
            $name = 'eid';
        }
        elsif ( $type eq 'notebook' ) {
            $name = 'nid';
        }

        if ( $format == 1 ) { # numbered
            $url .= join( ';',
                map { $name . ( $num++ ) . '=' . $_ } @{ $fields{$type} } );
        }
        elsif ( $format == 2 ) { # list
            $url .= $name . '=' . join( ',', @{ $fields{$type} } );
        }
        else {
            $url .= join( ';', map { $name . '=' . $_ } @{ $fields{$type} } );
        }
    }

    $url = $page_name . '.pl?' . $url if ($url);

    return $url;
}

sub get_contents {
    my %opts = @_;
    my $type = $opts{item_type};
    return unless $type;
    my $timestamp = $opts{timestamp};
    my $last_update = 0;
    #print STDERR "get_contents $type\n";
    
    #my $t1 = new Benchmark;
    my $items = [];

    if ( $type eq 'genome' ) {
        $items = get_genomes_for_user($DB->storage->dbh, $USER->id);
    }
    elsif ( $type eq 'experiment' ) {
        $items = get_experiments_for_user($DB->storage->dbh, $USER->id);
    }
    elsif ( $type eq 'notebook' ) {
        $items = get_lists_for_user($DB->storage->dbh, $USER->id);
    }
    elsif ( $type eq 'favorite' ) { # for favorites not owned by user (as those will be fetched in the above queries)
        my $favorites = CoGe::Core::Favorites->new(user => $USER);
        foreach ($favorites->get(notMine => 1)) {
            next if ($_->deleted);
            my $item = $_->to_hash();
            $item->{deleted}   = '0';
            $item->{favorite}  = '1';
            $item->{role_id}   = '2'; # fake an ownership role
            push @$items, $item;
        }
    }
    elsif ( $type eq 'metadata' ) {
		my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
		$template->param(
			METADATA => 1,
			EXPERIMENT_METADATA_STATS => get_stats('experiment', get_experiments_for_user($DB->storage->dbh, $USER->id)),
			GENOME_METADATA_STATS     => get_stats('genome', get_genomes_for_user($DB->storage->dbh, $USER->id)),
			NOTEBOOK_METADATA_STATS   => get_stats('list', get_lists_for_user($DB->storage->dbh, $USER->id))
		);
		return $template->output;
    }
    elsif ( $type eq 'group' ) {
        $items = get_groups_for_user($DB->storage->dbh, $USER->id);
    }
    elsif ( $type eq 'activity' ) {
        my $jobs = get_jobs($last_update);
        return summarize_jobs($jobs);
    }
    elsif ( $type eq 'analyses' ) {
        my $jobs = get_jobs($last_update);
        my $analyses = filter_jobs($jobs, ['synmap', 'cogeblast', 'gevo', 'synfind']);
        foreach (@$analyses) {
            push @$items, {
                id   => $_->{id},
                date => $_->{start_time},
                description => $_->{description},
                comment => $_->{comment},
                status => $_->{status},
                elapsed => $_->{elapsed},
                page => $_->{page},
                workflow_id => $_->{workflow_id},
                link => (is_uri($_->{link}) ? $_->{link} : undef),
                is_important => $_->{is_important}
            };
        }
    }
    elsif ( $type eq 'loads' ) {
        my $jobs = get_jobs($last_update);
        my $loads = filter_jobs($jobs, ['loadgenome', 'loadexperiment', 'loadannotation', 'loadbatch', 'genomeinfo', 'experimentview']);
        foreach (@$loads) {
            push @$items, {
                id   => $_->{id},
                date => $_->{start_time},
                description => $_->{description},
                status => $_->{status},
                elapsed => $_->{elapsed},
                page => $_->{page},
                workflow_id => $_->{workflow_id},
                link => (is_uri($_->{link}) ? $_->{link} : undef)
            };
        }
    }
    elsif ( $type eq 'graph' ) {
        # Generate activity graph
        my $user_id = $USER->id;
        my $job_list = 'cogeblast/synmap/gevo/synfind/loadgenome/loadexperiment/organismview/user';
        return qq{<iframe frameborder="0" width="100%" height="100%" scrolling="no" src="//genomevolution.org/blacktea/standalone/$user_id/$job_list"></iframe>} #FIXME: hardcoded server name
    }

#    print STDERR Dumper \@items, "\n";
    return encode_json($items);
}

sub get_stats {
	my $type = shift;
	my $items = shift;
	return '' if !$items;
    $items = [ grep { !$_->{'deleted'} } @$items ]; # if !$USER->is_admin();
	return '' if !scalar @$items;
	my $sql = 'select annotation_type.name,count(*) from ' . $type . '_annotation join annotation_type on annotation_type.annotation_type_id=' . $type . '_annotation.annotation_type_id where ' . $type . '_id in (' .
		join( ',', map { $_->{'id'} } @$items ) .
		') group by annotation_type.name';
	my $sth = $DB->storage->dbh->prepare($sql);
    $sth->execute();
    my ($name, $count);
    $sth->bind_columns(\$name, \$count);
    my @rows;
    while (my $row = $sth->fetch) {
        push @rows, [$name, $count];
    }
    my $num_rows = @rows;
    return $num_rows ? encode_json(\@rows) : 'null';
}

sub get_jobs {
    my ($last_update) = @_;
    
    # Get all jobs from the log table
    my @entries = $USER->logs(
        {   parent_id => { "!=" => undef},
            -or => [ parent_type => 7, parent_type => undef ], #FIXME magic number
            type => { '!=' => 0 },
            #time => { '>=' => $last_update } # for faster refresh
        },
        { order_by => { -desc => 'time' } }
    );
    
    # Add status info from job engine for currently active jobs
    my @workflow_ids = map { $_->parent_id } @entries;
    my $workflows = $JEX->find_workflows(\@workflow_ids);
    my %workflowsByID;
    foreach (@{$workflows}) {
        my($id, undef, $submitted, $completed, $status) = @{$_};
        my $diff = (defined $completed and defined $submitted) ? format_time_diff($completed - $submitted) : 0;
        $workflowsByID{$id} = {
            status => $status,
            elapsed => $diff,
        };
    }
    
    # Extract relevant fields
    my @results;
    foreach my $entry (@entries) {
        my $wid = $entry->parent_id;
        next unless defined $workflowsByID{$wid};
        
        push @results, {
            id => int($entry->id),
            page => $entry->page,
            description => $entry->description,
            comment => $entry->comment // '', #/
            link => $entry->link,
            start_time => $entry->time,
            is_important => $entry->is_important,
            workflow_id => $wid,
            elapsed => $workflowsByID{$wid}{elapsed},
            status => $workflowsByID{$wid}{status}
        };
    }

    # Add in old jobs from legacy jobs table
    foreach my $job ($USER->jobs({},
            {   time => { '>=' => $last_update },
                join => 'log',
                prefetch => 'log',
                order_by => { -desc => 'start_time' }
            }))
    {
        my $status = $job->status_description;
        $status =~ s/Running/Cancelled/; # Impossible to be running
        
        my $log = $job->log;
        push @results, {
            id => int($job->log_id),
            page => $job->page,
            description => $log->description,
            comment => $log->comment // '', #/
            link => $job->link,
            start_time => $job->start_time,
            is_important => $log->is_important,
            elapsed => $job->elapsed_time,
            status => $status
        };
    }
    
    # Filter out redundant results (kludge)
    my @filtered;
    foreach (reverse @results) {
        my $wid = $_->{workflow_id};
        next if (defined $wid and defined $workflowsByID{$wid}{seen});
        $workflowsByID{$wid}{seen}++ if (defined $wid);
        
        unshift @filtered, $_;
    }
    
    return \@filtered;
}

sub filter_jobs {
    my ($jobs, $page_names) = @_;
    
    my @results;
    foreach my $job (@$jobs) {
        next unless (grep { $job->{page} =~ /$_/i } @$page_names);
        push @results, $job;
    }
    
    return \@results;
}

sub summarize_jobs {
    my $jobs = shift;
    
    # Organize jobs by type and page
    my %page_counts;
    my %type_counts = ( loads => 0, analyses => 0 );
    foreach (@$jobs) {
        my $page = lc($_->{page});
        $page_counts{$page}++;
        if ( grep { $page =~ /$_/i } ('synmap', 'cogeblast', 'gevo', 'synfind') ) {
            $type_counts{analyses}++;
        }
        elsif ( grep { $page =~ /$_/i } ('loadgenome', 'loadexperiment', 'loadannotation', 'loadbatch') ) {
            $type_counts{loads}++;
        }
    }
    
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        ACTIVITY_SUMMARY => 1,
        NUM_ANALYSES => $type_counts{analyses},
        NUM_COGEBLAST => $page_counts{cogeblast} // 0, #/
        NUM_GEVO => $page_counts{gevo} // 0, #/
        NUM_SYNFIND => $page_counts{synfind} // 0, #/
        NUM_SYNMAP => $page_counts{synmap} // 0, #/
    );
    return $template->output;
}

sub format_job_status {
    my $status = shift;
    my $color;
    
    $status =~ s/terminated/cancelled/i;
    $status = lc($status);
    
    if ($status eq 'running')      { $color = 'yellowgreen' }
    elsif ($status eq 'completed') { $color = 'cornflowerblue' }
    elsif ($status eq 'scheduled') { $color = 'goldenrod' }
    else                           { $color = 'salmon' }
    
    return '<span style="padding-bottom:1px;padding-right:5px;padding-left:5px;border-radius:15px;color:white;background-color:' . $color . ';">' . ucfirst($status) . '</span>';
}

sub upload_image_file {
    return if ( $USER->user_name eq "public" );

    my %opts           = @_;
    my $image_filename = '' . $FORM->param('input_upload_file');
    my $fh             = $FORM->upload('input_upload_file');
    return if ( -s $fh > 2 * 1024 * 1024 ); # limit to 2MB

    #TODO delete old image

    # Create the image
    my $image;
    if ($fh) {
        $image = create_image(fh => $fh, filename => $image_filename, db => $DB);
        return 0 unless $image;

        # Link to user
        $USER->image_id( $image->id );
        $USER->update;
        return encode_json( { link => 'image.pl?id=' . $image->id } );
    }

    return;
}

sub upload_metadata {
    my %opts = @_;
    my $type = $opts{type};
    my @error_ids;
    my $file = $FORM->param('metadata_file');

    open my $fh, $FORM->tmpFileName($file);
    local $/ = undef;
	my $content = <$fh>;
	close $fh;

	my @lines = split /\r\n|\n|\r/, $content;
    my $line = shift @lines;
    chomp $line;
    my @headers = split(/\t/, $line);
    return "File doesn't seem to be tab delimited." if scalar @headers < 2;
    foreach (@lines) {
		chomp $_;
		next if !$_;
		my @row = split(/\t/, $_);
		my $annotations;
		for my $i (1..$#row) {
			next if !$row[$i];
			$annotations .= ';' if $annotations;
			$annotations .= $headers[$i] . '|' . $row[$i];
		}
		if ($annotations) {
			my @ids = split(/,/, $row[0]);
			foreach (@ids) {
			    my $can_annotate = 0;
			    if ($type eq 'Experiment') {
			        $can_annotate = $USER->is_owner_editor(experiment => $_);
			    }
			    elsif ($type eq 'Notebook') {
			        $can_annotate = $USER->is_owner_editor(list => $_);
			    }
			    else {
			        $can_annotate = $USER->is_owner_editor(dsg => $_);
			    }
			    if ($can_annotate) {
				    my $target = $DB->resultset($type eq 'Notebook' ? 'List' : $type)->find($_);
				    create_annotations(db => $DB, target => $target, annotations => $annotations, locked => 1);
			    }
			    else {
			        push @error_ids, $_;
			    }
			}
		}
	}
	if (@error_ids) {
	    return 'Permission denied for the following ' . $type . 's: ' . join(',', @error_ids);
	}
	return 'Metadata added';
}

sub search_notebooks
{    # FIXME this coded is dup'ed in CoGeBlast.pl and NotebookView.pl
    my %opts = @_;
    return if ( $USER->user_name eq 'public' );
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$search_term $timestamp\n";

    my @notebooks;
    my $num_results;
    my $group_str = join( ',', map { $_->id } $USER->groups );

    # Try to get all items if blank search term
    if ( !$search_term ) {
        my $sql = "locked=0"
          ;    # AND restricted=0 OR user_group_id IN ( $group_str ))"; # FIXME
        $num_results = $DB->resultset("List")->count_literal($sql);
        if ( $num_results < $MAX_SEARCH_RESULTS ) {
            foreach
              my $notebook ( $DB->resultset("List")->search_literal($sql) )
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
            $DB->resultset("List")->search_literal("locked=0 AND (name LIKE '$search_term' OR description LIKE '$search_term')")
        ) {
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
                html => "<option>$num_results matches, please refine your search.</option>"
            }
        );
    }

    # Build select items out of results
    my $html;
    foreach my $n ( sort notebookcmp @notebooks ) {
        $html .= "<option value='" . $n->id . "'>" . $n->info . "</option><br>\n";
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

    #print STDERR "add_items_to_notebook $nid $item_list\n";

    my $notebook = $DB->resultset('List')->find($nid);
    return unless $USER->has_access_to_list($notebook);

    foreach (@items) {
        my ( $item_id, $item_type ) = $_ =~ /(\d+)_(\w+)/;
        next unless ( $item_id and $item_type );
        next
          unless ( $item_type eq 'notebook'
            or $item_type eq 'genome'
            or $item_type eq 'experiment' );

        $item_type = $ITEM_TYPE{$item_type};
        #TODO check access permission on each item

        # print STDERR "add_item_to_notebook $item_id $item_type\n";

        my $conn = $DB->resultset('ListConnector')->find_or_create({
            parent_id  => $nid,
            child_id   => $item_id,
            child_type => $item_type
        });
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

    my $role = $DB->resultset('Role')->find($role_id);
    return unless $role;

    # Create the new group
    my $group = $DB->resultset('UserGroup')->create(
        {
            creator_user_id => $USER->id,
            name            => $name,
            description     => $desc,
            role_id         => $role->id
        }
    );
    return unless $group;

    # Set user as owner
    my $conn = $DB->resultset('UserConnector')->create(
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
		db          => $DB,
		user_id     => $USER->id,
		page        => $PAGE_TITLE,
		description => 'created group ' . $group->info_html,
		parent_id   => $group->id,
		parent_type => 6 #FIXME magic number
	);

    return 1;
}

sub create_new_notebook {
    my %opts    = @_;
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    #my $type_id = $opts{type_id}; # mdb removed 12/14/16 COGE-800
    return unless $name;
    my $item_list = $opts{item_list};    # optional
    return if ( $USER->user_name eq "public" );

    # Create the new list
    my $list = $DB->resultset('List')->create(
        {
            name         => $name,
            description  => $desc,
            #list_type_id => $type_id, # mdb removed 12/14/16 COGE-800
            creator_id   => $USER->id,
            restricted   => 1
        }
    );
    return unless $list;

    # Set user as owner
    my $conn = $DB->resultset('UserConnector')->create(
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
        db          => $DB,
        user_id     => $USER->id,
        page        => "$PAGE_TITLE",
        description => 'created notebook ' . $list->info_html,
        parent_id   => $list->id,
        parent_type => 1 #FIXME magic number
    );

    return 1;
}

# mdb removed 12/14/16 COGE-800
#sub get_notebook_types {
#    my $selected = shift;
#    my $html;
#    foreach my $type ( $DB->resultset('ListType')->all() ) {
#        next
#          if ( $type->name =~ /owner/i )
#          ;    # reserve this type for system-created lists
#        my $name =
#          $type
#          ->name;    # . ($type->description ? ": " . $type->description : '');
#        $html .=
#            '<option value="'
#          . $type->id . '" '
#          . (    $type->id eq $selected
#              || $type->name =~ /$selected/i ? 'selected' : '' )
#          . '>'
#          . $name
#          . '</option>';
#    }
#    return $html;
#}

sub toggle_star {
    my %opts   = @_;
    my $log_id = $opts{log_id};

    my $entry = $DB->resultset('Log')->find($log_id);
    return unless $entry;

    my $status = $entry->status;
    $entry->status( not $status );
    $entry->update();

    return not $status;
}

sub groupcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->name cmp $b->name;
}

sub usercmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->display_name cmp $b->display_name;
}
