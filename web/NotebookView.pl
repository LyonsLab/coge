#! /usr/bin/perl -w

use strict;
use CGI;
use JSON::XS;
use HTML::Template;
use Sort::Versions;
use List::Util qw(first);
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Metadata;
use CoGe::Core::Experiment qw(experimentcmp);
use CoGe::Core::Genome qw(genomecmp);
use CoGe::Core::Notebook qw(notebookcmp delete_notebook undelete_notebook);
use CoGe::Core::Storage qw(data_type);
use CoGeDBI qw(get_table_count);
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;
no warnings 'redefine';

use vars
  qw($P $PAGE_TITLE $PAGE_NAME $TEMPDIR $USER $BASEFILE $DB %FUNCTION $EMBED
  $FORM $TEMPDIR $TEMPURL $MAX_SEARCH_RESULTS $LINK $node_types);

$PAGE_TITLE = 'NotebookView';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $DB, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$TEMPDIR = $P->{TEMPDIR} . "$PAGE_TITLE/";
$TEMPURL = $P->{TEMPURL} . "$PAGE_TITLE/";

$MAX_SEARCH_RESULTS = 1000;
use constant MAX_TITLE_LENGTH => 150;

$node_types = $DB->node_types();

%FUNCTION = (
    get_list_info              => \&get_list_info,
    get_list_contents          => \&get_list_contents,
    edit_list_info             => \&edit_list_info,
    update_list_info           => \&update_list_info,
    make_list_public           => \&make_list_public,
    make_list_private          => \&make_list_private,
    add_list_items             => \&add_list_items,
    add_item_to_list           => \&add_item_to_list,
    remove_list_item           => \&remove_list_item,
    toggle_favorite            => \&toggle_favorite,
    get_annotations            => \&get_annotations,
    add_annotation             => \&add_annotation,
    get_annotation             => \&get_annotation,
    update_annotation          => \&update_annotation,
    remove_annotation          => \&remove_annotation,
    search_mystuff             => \&search_mystuff,
    search_genomes             => \&search_genomes,
    search_experiments         => \&search_experiments,
    search_features            => \&search_features,
    search_lists               => \&search_lists,
    search_users               => \&search_users,
    search_annotation_types    => \&search_annotation_types,
    get_annotation_type_groups => \&get_annotation_type_groups,
    delete_list                => \&delete_list,
    send_to_genomelist         => \&send_to_genomelist,
    send_to_experimentlist     => \&send_to_experimentlist,
    send_to_featlist           => \&send_to_featlist,
    send_to_blast              => \&send_to_blast,
    send_to_msa                => \&send_to_msa,
    send_to_gevo               => \&send_to_gevo,
    send_to_synfind            => \&send_to_synfind,
    send_to_featmap            => \&send_to_featmap,
    send_to_codeon             => \&send_to_codeon,
    send_to_fasta              => \&send_to_fasta,
    send_to_tsv                => \&send_to_tsv,
    send_to_xls                => \&send_to_xls,
    update_owner               => \&update_owner,
);

# debug for fileupload:
#print STDERR $ENV{'REQUEST_METHOD'} . "\n" . $FORM->url . "\n" . Dumper($FORM->Vars) . "\n";	# debug
#print "data begin\n" . $FORM->param('POSTDATA') . "\ndata end\n" if ($FORM->param('POSTDATA'));

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed') || 0;
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param(
            USER       => $USER->display_name || '',
            PAGE_TITLE => $PAGE_TITLE,
            TITLE      => "NotebookView",
            PAGE_LINK  => $LINK,
            SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
            HOME       => $P->{SERVER},
            HELP       => 'NotebookView',
            WIKI_URL   => $P->{WIKI_URL} || ''
        );
        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
        $template->param( ADMIN_ONLY => $USER->is_admin,
                          CAS_URL    => $P->{CAS_URL} || '',
                          COOKIE_NAME => $P->{COOKIE_NAME} || '' );
    }

    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $lid = $FORM->param('lid');
    $lid = $FORM->param('nid') unless $lid;                 # alias
    return "Must have valid notebook id\n" unless ($lid);
    my ($list) = $DB->resultset('List')->find($lid);
    return "<br>Notebook id$lid does not exist.<br>" unless ($list);
    return "Access denied\n" unless $USER->has_access_to_list($list);

    my $title = $list->info;
    $title = substr($title, 0, MAX_TITLE_LENGTH) . '...' if (length($title) > MAX_TITLE_LENGTH);

    my $favorites = CoGe::Core::Favorites->new(user => $USER);

    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        MAIN         => 1,
        EMBED        => $EMBED,
        PAGE_NAME    => $PAGE_TITLE . '.pl',
        NOTEBOOK_ID  => $lid,
        DEFAULT_TYPE => 'note',
        API_BASE_URL => $P->{SERVER} . 'api/v1/', #TODO move into config file or module
        USER         => $USER->user_name,
        NOTEBOOK_TITLE => $title,
        FAVORITED      => int($favorites->is_favorite($list)),
    );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( LIST_INFO => get_list_info( lid => $lid ) );
    $template->param( LIST_ANNOTATIONS => get_annotations( lid => $lid ) );
    $template->param( LIST_CONTENTS => get_list_contents( lid => $lid ) );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

sub get_list_info {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my ($list) = $DB->resultset('List')->find($lid);
    return unless $USER->has_access_to_list($list);

    my $html            = $list->annotation_pretty_print_html();
    my $user_can_edit   = $list->is_editable($USER);
    my $user_can_delete = $list->is_deletable($USER);

    $html .= qq{<div class="panel">};
    if ($user_can_edit) {
        $html .= qq{<span class='coge-button' style="margin-right:5px;" onClick="edit_list_info();">Edit Info</span>};

        if ( $list->restricted ) {
            $html .= qq{<span class='coge-button' style="margin-right:5px;" onClick="make_list_public();">Make Public</span>};
        }
        else {
            $html .= qq{<span class='coge-button' style="margin-right:5px;" onClick="make_list_private();">Make Private</span>};
        }
    }

    if ( $user_can_delete ) {
        $html .=
			qq{<span class='coge-button coge-button-danger' style="margin-right:5px;" onClick="delete_list();">} .
			($list->deleted ? 'Undelete' : 'Delete') . qq{</span>};
    }

    if ( $list->experiments( count => 1 ) ) {
        foreach my $gid (
            sort { $a <=> $b }
            map  { $_->genome_id } $list->experiments
          )
        {    # Pick a genome, any genome # TODO show user a list of genomes to choose from
            my $link = qq{window.open('GenomeView.pl?embed=$EMBED&gid=$gid&tracks=notebook$lid', '_self');};
            $html .= qq{<span class='coge-button' style="margin-right:5px;" onClick="$link">Browse</span>};
            last;
        }
    }
    
    $html .= qq{</div>};

    return $html;
}

# mdb removed 12/14/16 COGE-800
#sub get_list_types {
#    my $current_type_id = shift;
#
#    my @types;
#    foreach my $type ( $DB->resultset('ListType')->all() ) {
#        next if ( $type->name =~ /owner/i ); # reserve this type for system-created lists
#        my $name = $type->name . ( $type->description ? ": " . $type->description : '' );
#        my $selected = '';
#        $selected = 'selected="selected"' if ( $type->id == $current_type_id );
#        push @types,
#          { TID => $type->id, NAME => $name, TYPE_SELECTED => $selected };
#    }
#    return \@types;
#}

sub edit_list_info {
    my %opts = @_;
    my $lid  = $opts{lid};
    return 0 unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;

    my $desc = ( $list->description ? $list->description : '' );

    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        EDIT_LIST_INFO => 1,
        NAME           => $list->name,
        DESC           => $desc,
        #TYPE_LOOP      => get_list_types( $list->type->id ) # mdb removed 12/14/16 COGE-800
    );

    my %data;
    $data{name}   = $list->name;
    $data{desc}   = $desc;
    $data{output} = $template->output;

    return encode_json( \%data );
}

sub update_list_info {
    my %opts = @_;
    my $lid  = $opts{lid};
    return 0 unless $lid;
    my $name = $opts{name};
    return 0 unless $name;
    my $desc = $opts{desc};
    #my $type = $opts{type}; # mdb removed 12/14/16 COGE-800

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;

    $list->name($name);
    $list->description($desc) if $desc;
    #$list->list_type_id($type); # mdb removed 12/14/16 COGE-800
    $list->update;

    return 1;
}

sub make_list_public {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    #return unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;

    $list->restricted(0);
    $list->update;

    return 1;
}

sub make_list_private {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    #return unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;

    $list->restricted(1);
    $list->update;

    return 1;
}

sub toggle_favorite {
    my %opts   = @_;
    my $nid = $opts{nid};
    return unless $nid;
    return if ($USER->is_public); # must be logged in
    
    # Get genome
    my $notebook = $DB->resultset('List')->find($nid);
    return unless $notebook;
    
    # Toggle favorite
    my $favorites = CoGe::Core::Favorites->new( user => $USER );
    my $is_favorited = $favorites->toggle($notebook);
    
    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => ($is_favorited ? 'Favorited' : 'Unfavorited') . ' notebook ' . $notebook->info_html,
        parent_id   => $nid,
        parent_type => 1 #FIXME magic number
    );
    
    return $is_favorited;
}

sub linkify {
    my ( $link, $desc ) = @_;
    return "<span class='small link' onclick=\"window.open('$link')\">" . $desc . "</span>";
}

sub get_annotations {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless ($lid);
    my ($list) = $DB->resultset('List')->find($lid);
    return unless $USER->has_access_to_list($list);

    my $user_can_edit = $list->is_editable($USER);

    my %groups;
    my $num_annot = 0;
    foreach my $a ( $list->annotations ) {
        my $group = (
            defined $a->type->group
            ? $a->type->group->name . ':' . $a->type->name
            : $a->type->name
        );
        push @{ $groups{$group} }, $a;
        $num_annot++;
    }
    return unless ( $num_annot or $user_can_edit );

    my $html;
    if ($num_annot) {
        my $num_rows = 0;
        $html .= '<table id="list_annotation_table" class="dataTable compact hover stripe border-top border-bottom" style="max-width:400px;overflow:hidden;word-wrap:break-word;border-spacing:0;"><thead style="display:none"></thead><tbody>';
        foreach my $group ( sort keys %groups ) {
            my $first = 1;
            foreach my $a ( sort { $a->id <=> $b->id } @{ $groups{$group} } ) {
                $html .= "<tr style='vertical-align:top;' class='" . ($num_rows++ % 2 ? 'odd' : 'even') . "'>"
                  . (
                    $first-- > 0
                    ? "<th align='right' class='title5' rowspan='"
                      . @{ $groups{$group} }
                      . "' style='padding-right:10px;white-space:nowrap;font-weight:normal;background-color:white;'>$group:</th>"
                    : ''
                  );

                $html .= "<td>";
                my $image_link =
                  ( $a->image ? 'image.pl?id=' . $a->image->id : '' );
                my $image_info = (
                    $a->image
                    ? "<a href='$image_link' target='_blank' title='click for full-size image'><img height='40' width='40' src='$image_link' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='padding:1px;border:1px solid lightgray;vertical-align:text-top;'></a>"
                    : ''
                );
                $html .= $image_info if $image_info;
                $html .= "</td>";

                $html .= "<td class='data5'>" . $a->info . "</td>";
                $html .= "<td style='padding-left:5px;'>";
                $html .= linkify( $a->link, "Link" ) if $a->link;
                $html .= "</td>";
                if ($a->locked) {
                    $html .=
                        '<td style="padding-left:20px;white-space:nowrap;">'
                      . "<span onClick=\"alert('This item is locked and cannot be edited or removed.');\" class='link ui-icon ui-icon-locked'></span>"
                      . '</td>';
                }
                elsif ($user_can_edit) {
                    my $aid = $a->id;
                    $html .=
                        '<td style="padding-left:20px;white-space:nowrap;">'
                      . "<span onClick=\"edit_annotation_dialog($aid);\" class='link ui-icon ui-icon-gear'></span>"
                      . "<span onClick=\"\$(this).fadeOut(); remove_annotation($aid);\" class='link ui-icon ui-icon-trash'></span>"
                      . '</td>';
                }
                $html .= '</tr>';
            }
        }
        $html .= '</tbody></table>';
    }
    elsif ($user_can_edit) {
        $html .= '<table class="border-top border-bottom small padded note"><tr><td>There are no additional metadata items for this notebook.</tr></td></table>';
    }

    if ($user_can_edit) {
        $html .= qq{<div class="panel"><span onClick="add_annotation_dialog();" class='coge-button'>Add</span></div>};
    }

    return $html;
}

sub get_annotation {
    my %opts = @_;
    my $aid  = $opts{aid};
    return unless $aid;

    #TODO check user access here

    my $ea = $DB->resultset('ListAnnotation')->find($aid);
    return unless $ea;

    my $type       = '';
    my $type_group = '';
    if ( $ea->type ) {
        $type = $ea->type->name;
        $type_group = $ea->type->group->name if ( $ea->type->group );
    }
    return encode_json(
        {
            annotation => $ea->annotation,
            link       => $ea->link,
            type       => $type,
            type_group => $type_group
        }
    );
}

sub add_annotation {
    my %opts = @_;
    my $lid  = $opts{parent_id};
    return 0 unless $lid;
    my $type_group = $opts{type_group};
    my $type       = $opts{type};
    return 0 unless $type;
    my $annotation     = $opts{annotation};
    my $link           = $opts{link};
    my $image_filename = $opts{edit_annotation_image};
    my $fh             = $FORM->upload('edit_annotation_image');

   #return "Image file is too large (>10MB)" if (-s $fh > 10*1024*1024); # FIXME
   #print STDERR "add_annotation: $lid\n";

    if ($link) {
        $link =~ s/^\s+//;
        $link = 'http://' . $link if ( !$link =~ /^(\w+)\:\/\// );
    }

    my $group_rs;
    if ($type_group) {
        $group_rs =
          $DB->resultset('AnnotationTypeGroup')
          ->find( { name => $type_group } );

        # Create type group if it doesn't already exist
        if ( !$group_rs ) {
            $group_rs =
              $DB->resultset('AnnotationTypeGroup')
              ->create( { name => $type_group } );
        }
    }

    my $type_rs;
    $type_rs = $DB->resultset('AnnotationType')->find(
        {
            name                     => $type,
            annotation_type_group_id => ( $group_rs ? $group_rs->id : undef )
        }
    );

    # Create type if it doesn't already exist
    if ( !$type_rs ) {
        $type_rs = $DB->resultset('AnnotationType')->create(
            {
                name => $type,
                annotation_type_group_id =>
                  ( $group_rs ? $group_rs->id : undef )
            }
        );
    }

    # Create the image
    my $image;
    if ($fh) {
        $image = create_image(fh => $fh, filename => $image_filename, db => $DB);
        return 0 unless $image;
    }

    # Create the annotation
    my $la = $DB->resultset('ListAnnotation')->create(
        {
            list_id            => $lid,
            annotation         => $annotation,
            link               => $link,
            annotation_type_id => $type_rs->id,
            image_id           => ( $image ? $image->id : undef )
        }
    );
    return 0 unless $la;

    return 1;
}

sub update_annotation {
    my %opts = @_;
    my $aid  = $opts{aid};
    return unless $aid;
    my $type_group = $opts{type_group};
    my $type       = $opts{type};
    return 0 unless $type;
    my $annotation     = $opts{annotation};
    my $link           = $opts{link};
    my $image_filename = $opts{edit_annotation_image};
    my $fh             = $FORM->upload('edit_annotation_image');

    #TODO check user access here

    my $ea = $DB->resultset('ListAnnotation')->find($aid);
    return unless $ea;

    # Create the type and type group if not already present
    my $group_rs;
    if ($type_group) {
        $group_rs = $DB->resultset('AnnotationTypeGroup')->find_or_create( { name => $type_group } );
    }
    my $type_rs = $DB->resultset('AnnotationType')->find_or_create({
        name                     => $type,
        annotation_type_group_id => ( $group_rs ? $group_rs->id : undef )
    });

    # Create the image
    #TODO if image was changed delete previous image
    my $image;
    if ($fh) {
        $image = create_image(fh => $fh, filename => $image_filename, db => $DB);
        return 0 unless $image;
    }

    $ea->annotation($annotation);
    $ea->link($link);
    $ea->annotation_type_id( $type_rs->id );
    $ea->image_id( $image->id ) if ($image);
    $ea->update;

    return;
}

sub remove_annotation {
    my %opts = @_;
    my $lid  = $opts{lid};
    return "No notebook ID specified" unless $lid;
    my $laid = $opts{laid};
    return "No notebook annotation ID specified" unless $laid;

#return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $list = $DB->resultset('List')->find($lid);
    return 0 if ( $list->locked && !$USER->is_admin );

    my $la = $DB->resultset('ListAnnotation')->find( { list_annotation_id => $laid } );
    return 0 unless $la;
    $la->delete();

    return 1;
}

sub get_list_contents {
    my %opts = @_;
    my $lid  = $opts{lid};

    return "Must have valid notebook id\n" unless ($lid);

    my $list = $DB->resultset('List')->find($lid);

    return "Notebook id$lid does not exist.<br>" unless $list;

    return "Access denied\n" unless $USER->has_access_to_list($list);

    my @genomes = grep { $USER->has_access_to_genome($_); } $list->genomes;
    my @experiments = grep { $USER->has_access_to_experiment($_); } $list->experiments;
    my @features = $list->features;

    my @rows;
    foreach my $genome ( sort genomecmp @genomes ) {
        my $info = $genome->info;
        $info =~ s/'/\\'/g;
        $info =~ s/\n/<br>/g;
        push @rows, [ 0, $genome->id, $info, $genome->date eq '0000-00-00 00:00:00' ? undef : $genome->date ];
    }
    foreach my $experiment ( sort experimentcmp @experiments ) {
        my $info = $experiment->info;
        $info =~ s/'/\\'/g;
        $info =~ s/\n/<br>/g;
        push @rows, [ 1, $experiment->id, $info, $experiment->date eq '0000-00-00 00:00:00' ? undef : $experiment->date ];
    }
    foreach my $feature ( sort featurecmp @features ) {
        my $info = $feature->info;
        $info =~ s/'/\\'/g;
        $info =~ s/\n/<br>/g;
        push @rows, [ 2, $feature->id, $info ];
    }

    return encode_json({
        contents => \@rows,
        counts => [scalar @genomes, scalar @experiments, scalar @features],
        user_can_edit => $list->is_editable($USER) ? 'true' : 'false'
    });
}

sub add_list_items {
    my %opts = @_;
    my $lid  = $opts{lid};
    return 0 unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;

    my $desc = ( $list->description ? $list->description : '' );

    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( ADD_LIST_ITEMS => 1 );
    $template->param( NAME           => $list->name );
    $template->param( DESC           => $desc );

    # Setup dialog title and data fields
    my %data;
    $data{name}   = $list->name;
    $data{desc}   = $desc;
    $data{output} = $template->output;

    return encode_json( \%data );
}

sub add_item_to_list {
    my %opts = @_;
    my $lid  = $opts{lid};
    return 0 unless $lid;
    my $item_spec = $opts{item_spec};
    return 0 unless $item_spec;

    my ( $item_type, $item_id ) = split( /:/, $item_spec );
    #print STDERR "add_item_to_list: lid=$lid item_type=$item_type item_id=$item_id\n";

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;
    return 0 unless $list->is_editable($USER);

    my $lc =
      $DB->resultset('ListConnector')->find_or_create(
        {   parent_id => $lid,
            child_id => $item_id,
            child_type => $item_type
        });
    return 0 unless $lc;

    my $type_name = $DB->node_type_name($item_type);
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        page        => "NotebookView",
        description => "added $type_name id$item_id to notebook " . $list->info_html,
        link        => "NotebookView.pl?nid=$lid",
        parent_id   => $lid,
        parent_type => 1 #FIXME magic number
    );

    return 1;
}

sub remove_list_item {
    my %opts = @_;
    my $lid  = $opts{lid};
    return 0 unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;
    return 0 unless $list->is_editable($USER);

    my $item_type = $opts{item_type};
    my $item_id   = $opts{item_id};
    #print STDERR "remove_list_item: lid=$lid item_type=$item_type item_id=$item_id\n";

    my $lc =
      $DB->resultset('ListConnector')->find(
        {   parent_id => $lid,
            child_id => $item_id,
            child_type => $item_type
        });
    if ($lc) {
        $lc->delete();

        my $type_name = $DB->node_type_name($item_type);
        CoGe::Accessory::Web::log_history(
            db          => $DB,
            user_id     => $USER->id,
            page        => "NotebookView",
            description => "removed $type_name id$item_id from notebook " . $list->info_html,
            link        => "NotebookView.pl?nid=$lid",
            parent_id   => $lid,
            parent_type => 1 #FIXME magic number
        );
    }

    return 1;
}

sub search_mystuff {
    my %opts        = @_;
    my $lid         = $opts{lid};
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$lid $search_term\n";
    return 0 unless $lid;

    # Get items already in list
    my $list = $DB->resultset('List')->find($lid);

    my %exists;
    my $children = $list->children_by_type;
    foreach my $type ( keys %$children ) {
        foreach ( @{ $children->{$type} } ) {
            $exists{$type}{ $_->id }++;
        }
    }

    # Get my stuff
    my %mystuff;
    my $num_results = 0;

    my $type = $node_types->{experiment};
    foreach my $e ( $USER->experiments )
    {    #(sort experimentcmp $USER->experiments) {
        if ( !$search_term or $e->info =~ /$search_term/i ) {
            push @{ $mystuff{$type} }, $e;
            last if $num_results++ > $MAX_SEARCH_RESULTS;
        }
    }

    $type = $node_types->{genome};
    foreach my $g ( $USER->genomes ) {    #(sort genomecmp $USER->genomes) {
        if ( !$search_term or $g->info =~ /$search_term/i ) {
            push @{ $mystuff{$type} }, $g;
            last if $num_results++ > $MAX_SEARCH_RESULTS;
        }
    }

    # removed because Result/User.pm doesn't have features
    # $type = $node_types->{feature};
    # foreach my $f ( $USER->features ) {    #(sort featurecmp $USER->features) {
    #     if ( !$search_term or $f->info =~ /$search_term/i ) {
    #         push @{ $mystuff{$type} }, $f;
    #         last if $num_results++ > $MAX_SEARCH_RESULTS;
    #     }
    # }

# mdb: nested notebooks not supported
#    $type = $node_types->{list};
#    foreach my $l ( $USER->lists ) {       #(sort notebookcmp $USER->lists) {
#        next if ( $l->id == $list->id );    # can't add a list to itself!
#        next if ( $l->locked );             # exclude user's master list
#        if ( !$search_term or $l->info =~ /$search_term/i ) {
#            push @{ $mystuff{$type} }, $l;
#            last if $num_results++ > $MAX_SEARCH_RESULTS;
#        }
#    }

    # Limit number of results displayed
    if ( $num_results > $MAX_SEARCH_RESULTS ) {
        return encode_json(
            {
                timestamp => $timestamp,
                html => "<option disabled='disabled'>>$MAX_SEARCH_RESULTS results, please refine your search.</option>"
            }
        );
    }

    # Build select items out of results
    my $html;
    foreach my $type ( sort keys %mystuff ) {
        my ($type_name) = grep { $node_types->{$_} eq $type } keys %$node_types;
        $type_name = 'notebook' if ( $type_name eq 'list' );
        my $type_count = @{ $mystuff{$type} };
        $html .= "<optgroup label='" . ucfirst($type_name) . "s ($type_count)'>";
        foreach my $item ( @{ $mystuff{$type} } ) {
            my $disable =
              $exists{$type}{ $item->id } ? "disabled='disabled'" : '';
            my $item_spec = $type . ':' . $item->id;
            $html .=
                "<option $disable value='$item_spec'>"
              . $item->info
              . "</option>\n";
        }
        $html .= '</optgroup>';
    }
    $html = "<option disabled='disabled'>No matching items</option>"
      unless $html;

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub search_genomes {
    my %opts        = @_;
    my $lid         = $opts{lid};
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$lid $search_term $timestamp\n";
    return 0 unless $lid;

    # Get genomes already in list
    my $list = $DB->resultset('List')->find($lid);
    my %exists;
    map { $exists{ $_->id }++ } $list->genomes;

    my %unique;
    my $num_results;

    # Try to get all items if blank search term
    if ( !$search_term ) {

        # Get all genomes
        $num_results = $DB->resultset("Genome")->count;
        if ( $num_results < $MAX_SEARCH_RESULTS ) {
            my @genomes = $DB->resultset("Genome")->all;
            map {
                $unique{ $_->id } = $_
                  if ( !$_->deleted and $USER->has_access_to_genome($_) )
            } @genomes;
        }
    }

    # Perform search
    else {
        $search_term = '%' . $search_term . '%';

        # Get all matching organisms
        my @organisms = $DB->resultset("Organism")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $search_term ],
                [ 'description', $search_term ]
            ]
        );

        # Get all matching genomes
        my @genomes = $DB->resultset("Genome")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $search_term ],
                [ 'description', $search_term ]
            ]
        );

# Combine matching genomes with matching organism genomes, preventing duplicates
        map {
            $unique{ $_->id } = $_
              if ( !$_->deleted and $USER->has_access_to_genome($_) )
        } @genomes;
        foreach my $organism (@organisms) {
            map {
                $unique{ $_->id } = $_
                  if ( !$_->deleted and $USER->has_access_to_genome($_) )
            } $organism->genomes;
        }

        $num_results = keys %unique;
    }

    # Limit number of results displayed
    if ( $num_results > $MAX_SEARCH_RESULTS ) {
        return encode_json(
            {
                timestamp => $timestamp,
                html =>
"<option disabled='disabled'>$num_results results, please refine your search.</option>"
            }
        );
    }

    # Build select options out of results
    my $html;
    foreach my $g ( sort genomecmp values %unique ) {
        my $disable = $exists{ $g->id } ? "disabled='disabled'" : '';
        my $item_spec = $node_types->{genome} . ':' . $g->id;
        $html .=
          "<option $disable value='$item_spec'>" . $g->info . "</option><br>\n";
    }
    $html = "<option disabled='disabled'>No matching items</option>"
      unless $html;

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub search_experiments {
    my %opts        = @_;
    my $lid         = $opts{lid};
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$lid $search_term\n";
    return 0 unless $lid;

    # Get experiments already in list
    my $list = $DB->resultset('List')->find($lid);
    my %exists;
    map { $exists{ $_->id }++ } $list->experiments;

    my %unique;
    my $num_results;

    # Try to get all items if blank search term
    if ( !$search_term ) {

        # Get all experiments
        $num_results = $DB->resultset("Experiment")->count;
        if ( $num_results < $MAX_SEARCH_RESULTS ) {
            map { $unique{ $_->id } = $_ }
              $DB->resultset("Experiment")->search( { deleted => 0 } );
        }
    }

    # Perform search
    else {

        # Get user's private experiments
        foreach ( $USER->experiments( restricted => 1 ) ) {
            $unique{ $_->id } = $_
              if ( $_->name =~ /$search_term/i
                or $_->description =~ /$search_term/i );
        }

        # Get all public experiments
        $search_term = '%' . $search_term . '%';
        map { $unique{ $_->id } = $_ } $DB->resultset("Experiment")->search(
            \[
'restricted=? AND deleted=? AND (name LIKE ? OR description LIKE ?)',
                [ 'restricted',  0 ],
                [ 'deleted',     0 ],
                [ 'name',        $search_term ],
                [ 'description', $search_term ]
            ]
        );

        $num_results = keys %unique;
    }

    # Limit number of results displayed
    if ( $num_results > $MAX_SEARCH_RESULTS ) {
        return encode_json(
            {
                timestamp => $timestamp,
                html =>
"<option disabled='disabled'>$num_results results, please refine your search.</option>"
            }
        );
    }

    # Build select items out of results
    my $html;
    foreach my $exp ( sort experimentcmp values %unique ) {
        my $disable = $exists{ $exp->id } ? "disabled='disabled'" : '';
        my $item_spec = $node_types->{experiment} . ':' . $exp->id;
        $html .=
            "<option $disable value='$item_spec'>"
          . $exp->info
          . "</option><br>\n";
    }
    $html = "<option disabled='disabled'>No matching items</option>"
      unless $html;

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub search_features {
    my %opts        = @_;
    my $lid         = $opts{lid};
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$lid $search_term\n";
    return 0 unless $lid;

    # Get lists already in list
    my $list = $DB->resultset('List')->find($lid);
    my %exists;
    map { $exists{ $_->id }++ } $list->features;

    my @fnames;
    my $num_results;

    # Try to display all items if blank search term
    if ( !$search_term ) {

        # Get all features
        $num_results = get_table_count($DB->storage->dbh, 'feature_name');
        if ( $num_results < $MAX_SEARCH_RESULTS ) {
            @fnames = $DB->resultset("FeatureName")->all;
        }
    }

    # Perform search
    else {

        # Fulltext search copied from FeatView.pl
        push @fnames, $DB->resultset('FeatureName')->search( { name => $search_term } );
        unless (@fnames) {
            push @fnames,
              $DB->resultset('FeatureName')
              ->search_literal( 'MATCH(me.name) AGAINST (?)', $search_term );
        }
        $num_results = @fnames;
    }

    # Limit the number of results displayed
    if ( $num_results > $MAX_SEARCH_RESULTS ) {
        return encode_json(
            {
                timestamp => $timestamp,
                html =>
"<option disabled='disabled'>$num_results results, please refine your search.</option>"
            }
        );
    }

    # Build select items out of results
    my $html;
    my %seen;
    foreach my $f ( sort featurecmp @fnames ) {
        next if ( $seen{ $f->feature_id }++ );
        my $disable = $exists{ $f->feature_id } ? "disabled='disabled'" : '';
        my $item_spec = $node_types->{feature} . ':' . $f->feature_id;
        $html .=
            "<option $disable value='$item_spec'>"
          . $f->feature->info
          . "</option><br>\n";
    }
    $html = "<option disabled='disabled'>No matching items</option>"
      unless $html;

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub search_lists
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
            $DB->resultset("List")->search_literal(
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
    foreach my $n ( sort notebookcmp @notebooks ) {
        my $item_spec = 1 . ':' . $n->id;    #FIXME magic number for item_type
        $html .= "<option value='$item_spec'>" . $n->info . "</option><br>\n";
    }
    $html = "<option disabled='disabled'>No matches</option>" unless $html;

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub search_users {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @users = $DB->resultset("User")->search(
        \[
            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
            [ 'user_name',  $search_term ],
            [ 'first_name', $search_term ],
            [ 'last_name',  $search_term ]
        ]
    );

    # Limit number of results displayed
    # if (@users > $MAX_SEARCH_RESULTS) {
    #   return encode_json({timestamp => $timestamp, items => undef});
    # }

    return encode_json(
        {
            timestamp => $timestamp,
            items     => [ sort map { $_->user_name } @users ]
        }
    );
}

sub search_annotation_types {
    my %opts        = @_;
    my $type_group  = $opts{type_group};
    my $search_term = $opts{search_term};

    #	print STDERR "search_annotation_types: $search_term $type_group\n";
    return '' unless $search_term;

    $search_term = '%' . $search_term . '%';

    my $group;
    if ($type_group) {
        $group =
          $DB->resultset('AnnotationTypeGroup')
          ->find( { name => $type_group } );
    }

    my @types;
    if ($group) {

        #		print STDERR "type_group=$type_group " . $group->id . "\n";
        @types = $DB->resultset("AnnotationType")->search(
            \[
'annotation_type_group_id = ? AND (name LIKE ? OR description LIKE ?)',
                [ 'annotation_type_group_id', $group->id ],
                [ 'name',                     $search_term ],
                [ 'description',              $search_term ]
            ]
        );
    }
    else {
        @types = $DB->resultset("AnnotationType")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $search_term ],
                [ 'description', $search_term ]
            ]
        );
    }

    my %unique;
    map { $unique{ $_->name }++ } @types;
    return encode_json( [ sort keys %unique ] );
}

sub get_annotation_type_groups {

    #my %opts = @_;
    my %unique;

    my $rs = $DB->resultset('AnnotationTypeGroup');
    while ( my $atg = $rs->next ) {
        $unique{ $atg->name }++;
    }

    return encode_json( [ sort keys %unique ] );
}

sub featurecmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->name cmp $b->name;
}

sub delete_list {
    my %opts = @_;
    my $lid  = $opts{lid};
    return 0 unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return 0 unless $list;

    my $success;
    if ($list->deleted) {
        $success = undelete_notebook(db => $DB, user => $USER, notebook_id => $lid, page => $PAGE_TITLE);
    }
    else {
        $success = delete_notebook(db => $DB, user => $USER, notebook_id => $lid, page => $PAGE_TITLE);
    }
   
    return $success;
}

sub send_to_blast {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    #my $list = $DB->resultset('List')->find($lid);
    #return unless $list;

    #my $accn_list = join(',', map { $_->id } $list->genomes);
    #my $url = "CoGeBlast.pl?dsgid=$accn_list";
    my $url = "CoGeBlast.pl?lid=$lid";
    return encode_json( { url => $url } );
}

sub send_to_genomelist {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;
    my $url = "GenomeList.pl?lid=$lid";
    return encode_json( { url => $url } );
}

sub send_to_experimentlist {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;
    my $url = "ExperimentList.pl?lid=$lid";
    return encode_json( { url => $url } );
}

sub send_to_featlist {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;
    my $url = "FeatList.pl?lid=$lid";
    return encode_json( { url => $url } );
}

sub send_to_msa {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $accn_list = join( '&', map { 'fid=' . $_->id } $list->features );
    my $url = "CoGeAlign.pl?$accn_list";
    return encode_json( { url => $url } );
}

sub send_to_gevo {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $count = 1;
    my $accn_list =
      join( '&', map { 'fid' . $count++ . '=' . $_->id } $list->features );

    my $url = "GEvo.pl?$accn_list";
    $count--;
    return encode_json(
        {
            alert =>
"You have exceeded the number of features you can send to GEvo (20 max)."
        }
    ) if $count > 20;
    $url .= "&num_seqs=$count";
    return encode_json( { url => $url } );
}

sub send_to_synfind {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $accn_list = join( ',', map { $_->id } $list->genomes );
    my $url       = "SynFind.pl?dsgid=$accn_list";
    my $first     = first { $_->id } $list->features;
    $url .= ";fid=" . $first->id if ($first);
    return encode_json( { url => $url } );
}

sub send_to_featmap {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $accn_list = join( '&', map { 'fid=' . $_->id } $list->features );
    my $url = "FeatMap.pl?$accn_list";
    return encode_json( { url => $url } );
}

sub send_to_codeon {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $accn_list = join( '::', map { $_->id } $list->features );
    my $url = "CodeOn.pl?fid=$accn_list";
    return encode_json( { url => $url } );
}

sub send_to_fasta {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $cogeweb =
      CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $basename = $cogeweb->basefilename;
    my $file     = $TEMPDIR . "$basename.faa";
    open( OUT, ">$file" );
    foreach my $g ( sort genomecmp $list->genomes ) {
        print OUT $g->fasta;
    }
    close OUT;
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
}

sub send_to_tsv {
    my %opts = @_;
    my $lid  = $opts{lid};
    return unless $lid;

    my $list = $DB->resultset('List')->find($lid);
    return unless $list;

    my $cogeweb =
      CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $basename = $cogeweb->basefilename;
    my $file     = "$TEMPDIR/$basename.tsv";

    #	print STDERR $P->{SERVER} . " $TEMPURL $file\n";
    open( OUT, ">$file" );

    my @genomes = $list->genomes;
    if (scalar @genomes) {
        print OUT join( "\t",
            "CoGe Genome ID",
            "Name",
            "Description",
            "Source",
            "Provenance",
            "Sequence Type",
            "Chr Count",
            "Length (bp)",
            "Percent GC",
            "Percent AT",
            "Percent N|X",
            "OrganismView Link" ),
            "\n";
        foreach my $g ( sort genomecmp @genomes ) {
            my $name = $g->name ? $g->name : $g->organism->name;
            my $desc = $g->description ? $g->description : $g->organism->description;
            my ($ds_source) = $g->source;
            my $source = $ds_source->name;
            my $provenance = join( "||", map { $_->name } $g->datasets );
            my $chr_count  = $g->chromosome_count;
            my $length     = $g->length;
            my $type       = $g->type->name;
            warn 1;
            my ( $gc, $at, $n ) = $g->percent_gc();
            warn 2;
            $at *= 100;
            $gc *= 100;
            print OUT join( "\t",
                $g->id,      $name,
                $desc,       $source,
                $provenance, $type,
                $chr_count,  $length,
                $gc,         $at,
                $n,          $P->{SERVER} . 'OrganismView.pl?dsgid=' . $g->id ),
                "\n";
        }
    }

    my @experiments = $list->experiments;
    if (scalar @experiments) {
        print OUT join( "\t",
            "CoGe Experiment ID",
            "Name",
            "Description",
            "Data Type",
            "Genome ID",
            "Source",
            "Version" ),
            "\n";
        foreach my $e ( sort experimentcmp @experiments ) {
            print OUT join( "\t",
                $e->id,
                $e->name,
                $e->description,
                data_type($e->data_type),
                $e->genome_id,
                $e->source->name,
                $e->version ),
                "\n";
        }
    }

    close OUT;
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return encode_json( { url => $file } );
}

sub send_to_xls {
    my %args      = @_;
    my $accn_list = $args{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $cogeweb =
      CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $basename = $cogeweb->basefilename;
    my $file     = "$TEMPDIR/Excel_$basename.xls";
    my $workbook = Spreadsheet::WriteExcel->new($file);
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i         = 1;

    $worksheet->write( 0, 0,  "Name" );
    $worksheet->write( 0, 1,  "Description" );
    $worksheet->write( 0, 2,  "Source" );
    $worksheet->write( 0, 3,  "Provenance" );
    $worksheet->write( 0, 4,  "Sequence Type" );
    $worksheet->write( 0, 5,  "Chr Count" );
    $worksheet->write( 0, 6,  "Length (bp)" );
    $worksheet->write( 0, 7,  "Percent GC" );
    $worksheet->write( 0, 8,  "Percent AT" );
    $worksheet->write( 0, 9,  "Percent N|X" );
    $worksheet->write( 0, 10, "OrganismView Link" );

    foreach my $dsgid ( split( /,/, $accn_list ) ) {
        my ($dsg) = $DB->resultset("Genome")->find($dsgid);
        next unless $dsg;

        my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
        my $desc =
          $dsg->description ? $dsg->description : $dsg->organism->description;
        my ($ds_source) = $dsg->source;
        my $source = $ds_source->name;
        my $provenance = join( " ", map { $_->name } $dsg->datasets );
        my $length     = $dsg->length;
        my $chr_count  = $dsg->chromosome_count;
        my $type       = $dsg->type->name;
        my ( $gc, $at, $n ) = $dsg->percent_gc();
        $at *= 100;
        $gc *= 100;

        #my ($wgc, $wat) = $dsg->wobble_content;
        #$wat*=100;
        #$wgc*=100;

        $worksheet->write( $i, 0, $name );
        $worksheet->write( $i, 1, $desc );
        $worksheet->write( $i, 2, $source );
        $worksheet->write( $i, 3, $provenance );
        $worksheet->write( $i, 4, $type );
        $worksheet->write( $i, 5, $chr_count );
        $worksheet->write( $i, 6, $length );
        $worksheet->write( $i, 7, $gc . '%' );
        $worksheet->write( $i, 8, $at . '%' );
        $worksheet->write( $i, 9, $n . '%' );
        $worksheet->write( $i, 10,
            $P->{SERVER} . 'OrganismView.pl?dsgid=' . $dsgid );

        $i++;
    }
    $workbook->close() or die "Error closing file: $!";
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
}

sub update_owner {
    my %opts      = @_;
    my $nid       = $opts{nid};
    my $user_name = $opts{user_name};
    return unless $nid and $user_name;

    # Admin-only function
    return unless $USER->is_admin;

    # Get user to become owner
    my $user = $DB->resultset('User')->find( { user_name => $user_name } );
    unless ($user) {
        return "error finding user '$user_name'\n";
    }

    # Get notebook
    my $notebook = $DB->resultset('List')->find($nid);
    unless ($notebook) {
        return "error finding notebook id$nid\n";
    }    

    # Make user owner of notebook and its contents
    foreach my $item ($notebook, $notebook->children) {
        print STDERR 'Reassigning ', $item->id, ' ', $item->item_type, "\n";
        
        # Remove current owner (should only be one, but loop just in case)
        foreach my $conn (
            $DB->resultset('UserConnector')->search(
                {
                    parent_type => $node_types->{user},
                    child_id    => $item->id,
                    child_type  => $item->item_type,
                    role_id     => 2                        # FIXME hardcoded
                }
            ))
        {
            $conn->delete;
        }
        
        # Remove existing user connection (should only be one, but loop just in case)
        foreach my $conn (
            $DB->resultset('UserConnector')->search(
                {
                    parent_id   => $user->id,
                    parent_type => $node_types->{user},
                    child_id    => $item->id,
                    child_type  => $item->item_type,
                    role_id     => {'!=' => 2}           # FIXME hardcoded
                }
            ))
        {
            $conn->delete;
        }
        
        # Make given user the owner
        my $conn = $DB->resultset('UserConnector')->find_or_create(
            {
                parent_id   => $user->id,
                parent_type => $node_types->{user},
                child_id    => $item->id,
                child_type  => $item->item_type,
                role_id     => 2                        # FIXME hardcoded
            }
        );
        unless ($conn) {
            return "error creating user connector\n";
        }
    }

    return;
}