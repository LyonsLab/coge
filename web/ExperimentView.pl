#! /usr/bin/perl -w

use strict;

use CGI;
use HTML::Template;
use JSON::XS;
use Spreadsheet::WriteExcel;
use File::Basename;
use File::Path;
use File::Slurp;
use File::Spec::Functions qw(catdir catfile);
use Sort::Versions;
use Switch;
use Data::Dumper;

use CoGe::Accessory::IRODS qw(irods_get_base_path irods_imeta_add irods_iput);
use CoGe::Accessory::Utils;
use CoGe::Accessory::Web;
use CoGe::Core::Experiment qw(get_irods_metadata);
use CoGe::Core::Favorites;
use CoGe::Core::Metadata;
use CoGe::Core::Storage;

use vars qw(
    $P $PAGE_TITLE $USER $LINK $DB $FORM $EMBED %FUNCTION $ERROR $WORKFLOW_ID $LOAD_ID $TEMPDIR
);

$PAGE_TITLE = "ExperimentView";
$ERROR = encode_json( { error => 1 } );

use constant MAX_TITLE_LENGTH => 150;

$FORM = new CGI;
( $DB, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE,
);

$WORKFLOW_ID = $FORM->Vars->{'wid'} || $FORM->Vars->{'job_id'}; # wid is new name, job_id is legacy name
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = $P->{SECTEMPDIR} . $PAGE_TITLE . '/' . $USER->name . '/' . $LOAD_ID . '/';

%FUNCTION = (
    get_experiment_info        => \&get_experiment_info,
    edit_experiment_info       => \&edit_experiment_info,
    update_experiment_info     => \&update_experiment_info,
    get_sources                => \&get_sources,
    toggle_favorite            => \&toggle_favorite,
    make_experiment_public     => \&make_experiment_public,
    make_experiment_private    => \&make_experiment_private,
    add_tag_to_experiment      => \&add_tag_to_experiment,
    get_experiment_tags        => \&get_experiment_tags,
    remove_experiment_tag      => \&remove_experiment_tag,
    get_annotations            => \&get_annotations,
    add_annotation             => \&add_annotation,
    update_annotation          => \&update_annotation,
    remove_annotation          => \&remove_annotation,
    search_annotation_types    => \&search_annotation_types,
    check_login                => \&check_login,
    export_experiment_irods    => \&export_experiment_irods,
    get_file_urls              => \&get_file_urls,
    get_progress_log           => \&get_progress_log,
    get_load_log               => \&get_load_log,
    send_error_report          => \&send_error_report
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub edit_experiment_info {
    my %opts = @_;
    my $eid  = $opts{eid};
    return 0 unless $eid;

    my $exp = $DB->resultset('Experiment')->find($eid);
    my $desc = ( $exp->description ? $exp->description : '' );


    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );

    $template->param(
        EDIT_EXPERIMENT_INFO => 1,
        EXPERIMENT_ID        => $eid,
        NAME                 => $exp->name,
        DESC                 => $desc,
        SOURCE               => $exp->source->name,
        SOURCE_ID            => $exp->source->id,
        LINK                 => $exp->link,
        VERSION              => $exp->version,
    );

    my %data;
    $data{name}   = $exp->name;
    $data{desc}   = $desc;
    $data{output} = $template->output;

    return encode_json( \%data );
}

sub update_experiment_info {
    my %opts = @_;
    my $eid  = $opts{eid};
    return 0 unless $eid;
    my $name = $opts{name};
    return 0 unless $name;
    my $desc      = $opts{desc};
    my $source_id = $opts{source_id};
    my $link      = $opts{link};
    my $version   = $opts{version};

    my $exp = $DB->resultset('Experiment')->find($eid);
    $exp->name($name);
    $exp->description($desc) if $desc;
    $exp->version($version);
    $exp->link($link);
    $exp->data_source_id($source_id);
    $exp->update;

    return 1;
}

sub get_sources {
    #my %opts = @_;

    my @sources;
    foreach ( $DB->resultset('DataSource')->all() ) {
        push @sources, { 'label' => $_->name, 'value' => $_->id };
    }

    return encode_json( \@sources );
}

sub make_experiment_public {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No EID specified" unless $eid;
    #return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $exp = $DB->resultset('Experiment')->find($eid);
    $exp->restricted(0);
    $exp->update;

    return 1;
}

sub make_experiment_private {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No EID specified" unless $eid;
    #return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $exp = $DB->resultset('Experiment')->find($eid);
    $exp->restricted(1);
    $exp->update;

    return 1;
}

sub add_tag_to_experiment {
    my %opts = @_;
    my $eid  = $opts{eid};
    return 0 unless $eid;
    my $name = $opts{name};
    return 0 unless $name;

    my $type = $DB->resultset('ExperimentType')->find( { name => $name } );

    if ($type) {
        # If type exists, check if already assigned to this experiment
        foreach ( $type->experiment_type_connectors ) {
            return 1 if ( $_->experiment_id == $eid );
        }
    }
    else {
        # Create type if it doesn't already exist
        $type =
          $DB->resultset('ExperimentType')
          ->create( { name => $name } );
    }

    # Create connection
    $DB->resultset('ExperimentTypeConnector')->create(
        {
            experiment_id      => $eid,
            experiment_type_id => $type->id
        }
    );

    return 1;
}

#FIXME: Types should be more generic and be referred to as TAGS
sub get_experiment_tags {
    my $tags = $DB->storage->dbh->selectall_arrayref('SELECT DISTINCT name FROM experiment_type ORDER BY name');
    my @a = map { $_->[0] } @$tags;
    return encode_json(\@a);
}

sub linkify {
    my ( $link, $desc ) = @_;
    return "<span class='link' onclick=\"window.open('$link')\">$desc</span>";
}

sub remove_experiment_tag {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No experiment ID specified" unless $eid;
    my $etid = $opts{etid};
    return "No experiment type ID specified" unless $etid;
    #return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $etc =
      $DB->resultset('ExperimentTypeConnector')
      ->find( { experiment_id => $eid, experiment_type_id => $etid } );
    $etc->delete();

    return 1;
}

sub toggle_favorite {
    my %opts   = @_;
    my $eid = $opts{eid};
    return unless $eid;
    return if ($USER->is_public); # must be logged in
    
    # Get genome
    my $experiment = $DB->resultset('Experiment')->find($eid);
    return unless $experiment;
    
    # Toggle favorite
    my $favorites = CoGe::Core::Favorites->new( user => $USER );
    my $is_favorited = $favorites->toggle($experiment);
    
    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => ($is_favorited ? 'Favorited' : 'Unfavorited') . ' experiment ' . $experiment->info_html,
        parent_id   => $eid,
        parent_type => 3 #FIXME magic number
    );
    
    return $is_favorited;
}

sub get_annotations {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "Must have valid experiment id\n" unless ($eid);

    my $exp = $DB->resultset('Experiment')->find($eid);
    return "Access denied\n" unless $USER->has_access_to_experiment($exp);

    my $user_can_edit = $exp->is_editable($USER);

    # Categorize annotations based on type group and type
    my %groups;
    my $num_annot = 0;
    foreach my $a ( $exp->annotations ) {
        my $group = ( $a->type->group ? $a->type->group->name : '');
        my $type = $a->type->name;
        push @{ $groups{$group}{$type} }, $a if (defined $group and defined $type);
        $num_annot++;
    }

    # Build annotation table
    my $html;
    if ($num_annot) {
        $html .= '<table id="experiment_annotation_table" class="border-top border-bottom small" style="max-width:800px;overflow:hidden;word-wrap:break-word;border-spacing:0;"><thead style="display:none"></thead><tbody>';
        foreach my $group ( sort keys %groups ) { # groups
            my $first_group = 1;
            foreach my $type ( sort keys %{ $groups{$group} } ) { # types
                my $first_type = 1;
                foreach my $a ( sort { $a->id <=> $b->id } @{ $groups{$group}{$type} } ) { # annotations
                    my $header = ($group and $first_group-- > 0 ? "<b>$group</b> " : '') . ($first_type-- > 0 ? $type : '');
                    $html .= "<tr style='vertical-align:top;'>";
                    $html .= "<th class='title4' style='padding-right:10px;white-space:nowrap;font-weight:normal;background-color:white;text-align:left;'>$header</th><td class='data4'>";
                    if ($a->image) {
                        my $image_link = ( $a->image ? 'image.pl?id=' . $a->image->id : '' );
                        $html .= "<a href='$image_link' target='_blank' title='click for full-size image'><img height='40' width='40' src='$image_link' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='float:left;padding:1px;border:1px solid lightgray margin-right:5px;'></a>";
                    }
                    elsif ($a->bisque_id) {
                        $html .= "<a href='http://bisque.iplantc.org/client_service/view?resource=http://bisque.iplantc.org/data_service/";
                        $html .= $a->bisque_id;
                        $html .= "6' target='_blank' title='click to view in BisQue'><img src='http://bisque.iplantc.org/image_service/";
                        $html .= $a->bisque_id;
                        $html .= "6?thumbnail=200,200' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='float:left;padding:1px;border:1px solid lightgray;width:42px;margin-right:5px;'></a>";
                    }
                    $html .= $a->info;
                    $html .= '</td><td style="padding-left:5px;">';
                    $html .= linkify( $a->link, 'Link' ) if $a->link;
                    $html .= '</td>';
                    if ($user_can_edit && !$a->locked) {
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
        }
        $html .= '</tbody></table>';
    }
    elsif ($user_can_edit) {
        $html .= '<table class="border-top border-bottom small padded note"><tr><td>There are no additional metadata items for this experiment.</tr></td></table>';
    }

    if ($user_can_edit) {
        $html .= qq{<div class="padded"><span onClick="add_annotation_dialog();" class='coge-button'>Add</span></div>};
    }

    return $html;
}

sub add_annotation {
    my %opts = @_;
    my $fh = $FORM->upload('edit_annotation_image');
 
    return CoGe::Core::Metadata::create_annotation(
        db => $DB,
        fh => $fh,
        filename => $opts{edit_annotation_image},
        group_name => $opts{type_group},
        image_tmp_file => $FORM->tmpFileName($opts{edit_annotation_image}),
        link => $opts{link},
        locked => 0,
        target_id => $opts{parent_id},
        target_type => 'experiment',
        text => $opts{annotation},
        type_name => $opts{type},
        user => $USER
    ) ? 1 : 0;
}

sub update_annotation {
    my %opts = @_;
    my $fh = $FORM->upload('edit_annotation_image');
    CoGe::Core::Metadata::update_annotation(
        annotation_id => $opts{aid},
        db => $DB,
        fh => $fh,
        filename => $opts{edit_annotation_image},
        group_name => $opts{type_group},
        image_tmp_file => $FORM->tmpFileName($opts{edit_annotation_image}),
        link => $opts{link},
        target_type => 'experiment',
        text => $opts{annotation},
        type_name => $opts{type},
        user => $USER
    );
    return 1;
}

#XXX: Move to a module
sub check_login {
    #print STDERR $USER->user_name . ' ' . int($USER->is_public) . "\n";
    return ($USER && !$USER->is_public);
}

sub export_experiment_irods { #TODO migrate to API
    my %opts = @_;
    my $eid = $opts{eid};

    my $experiment = $DB->resultset('Experiment')->find($eid);
    return $ERROR unless $USER->has_access_to_experiment($experiment);

    my ($statusCode, $file) = generate_export($experiment);
    unless($statusCode) {
        my $dir = irods_get_base_path($USER->user_name);
        my $dest = File::Spec->catdir($dir, basename($file));

        my $md = get_irods_metadata($experiment);

        irods_iput($file, $dest);
        irods_imeta_add($dest, $md);
    }

    return basename($file);
}

sub generate_export { #TODO migrate to API
    my $experiment = shift;
    my $eid = $experiment->id;

    my $exp_name = sanitize_name($experiment->name);
       $exp_name = $eid unless $exp_name;

    my $filename = "experiment_$exp_name.tar.gz";

    my $conf = $P->{_CONFIG_PATH};
    my $script = File::Spec->catdir($P->{SCRIPTDIR}, "export_experiment_or_genome.pl");
    my $workdir = get_experiment_cache_path($eid);
    my $resdir = $P->{RESOURCEDIR};

    my $cmd = "$script -id $eid -type 'experiment' -config $conf -dir $workdir -output $filename";

    return (execute($cmd), File::Spec->catdir($workdir, $filename));
}

sub get_file_urls {
    my %opts = @_;
    my $eid = $opts{eid};
    return 0 unless $eid;

    my $experiment = $DB->resultset('Experiment')->find($eid);
    return 0 unless $USER->has_access_to_experiment($experiment);

    my ($statusCode, $file) = generate_export($experiment);

    unless($statusCode) {
        my $url = download_url_for(eid => $eid, file => $file);
        return encode_json({ filename => basename($file), url => $url });
    };

    return encode_json({ error => 1 });
}

sub remove_annotation {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No experiment ID specified" unless $eid;
    my $eaid = $opts{eaid};
    return "No experiment annotation ID specified" unless $eaid;
    #return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $ea = $DB->resultset('ExperimentAnnotation')->find( { experiment_annotation_id => $eaid } );
    $ea->delete();

    return 1;
}

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed') || 0;
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new(filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param(
            PAGE_TITLE => $PAGE_TITLE,
            TITLE      => 'ExperimentView',
            PAGE_LINK  => $LINK,
            SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
            HOME       => $P->{SERVER},
            HELP       => 'ExperimentView',
            WIKI_URL   => $P->{WIKI_URL} || '',
            CAS_URL    => $P->{CAS_URL} || '',
            COOKIE_NAME => $P->{COOKIE_NAME} || ''
        );
    	$template->param( USER     => $USER->display_name || '' );
        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
        $template->param( ADMIN_ONLY => $USER->is_admin );
    }

    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $eid = $FORM->param('eid');
    return "Need a valid experiment id\n" unless $eid;

    my $exp = $DB->resultset('Experiment')->find($eid);
    return "Experiment not found" unless $exp;
    return "Access denied" unless $USER->has_access_to_experiment($exp);

    my $gid = $exp->genome_id;
    
    my $popgenUrl;
    if ($exp->data_type == $DATA_TYPE_POLY && is_popgen_finished($eid)) {
        $popgenUrl = "PopGen.pl?eid=$eid";
    }

    my $title = $exp->info;
    $title = substr($title, 0, MAX_TITLE_LENGTH) . '...' if (length($title) > MAX_TITLE_LENGTH);

    my $favorites = CoGe::Core::Favorites->new(user => $USER);

    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        MAIN              => 1,
        EMBED             => $EMBED,
        PAGE_NAME         => $PAGE_TITLE . '.pl',
        USER_NAME         => $USER->name,
        EXPERIMENT_ID     => $eid,
        EXPERIMENT_TITLE  => $title,
        FAVORITED         => int($favorites->is_favorite($exp)),
        DEFAULT_TYPE      => 'note',
        ITEMS             => commify($exp->row_count),
        FILE_SIZE         => commify(directory_size(get_experiment_path($exp->id))),
        IRODS_HOME        => irods_get_base_path($USER->user_name),
        WORKFLOW_ID       => $WORKFLOW_ID,
        STATUS_URL        => 'jex/status/',
        ALIGNMENT_TYPE    => ($exp->data_type == $DATA_TYPE_ALIGN),
        POLYMORPHISM_TYPE => ($exp->data_type == $DATA_TYPE_POLY),
        POPGEN_RESULT_URL => $popgenUrl,
        PUBLIC            => $USER->user_name eq "public" ? 1 : 0,
        ADMIN_AREA        => $USER->is_admin,
        API_BASE_URL      => $P->{SERVER} . 'api/v1/', #TODO move into config file or module
    );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( EXPERIMENT_INFO => get_experiment_info( eid => $eid ) || undef );
    $template->param( EXPERIMENT_ANNOTATIONS => get_annotations( eid => $eid ) || undef );

    return $template->output;
}


sub _get_experiment_info {
    my %opts  = @_;
    my $eid   = $opts{eid};
    my ($exp) = $DB->resultset('Experiment')->find($eid);
    return "Access denied\n" unless $USER->has_access_to_experiment($exp);
    return "Unable to find an entry for $eid" unless $exp;

    my $allow_edit = $USER->is_admin || $USER->is_owner_editor( experiment => $eid );
    my $gid = $exp->genome->id;

    my $tags;
    foreach my $tag ( $exp->tags ) {
       $tags .= '<span class="coge-tag">';
       $tags .= $tag->name;
       if ($allow_edit) {
           $tags .=
               "<span onClick=\"remove_experiment_tag({eid: '"
             . $exp->id
             . "', etid: '"
             . $tag->id
             . "'});\" class=\"link ui-icon ui-icon-close\"></span>";
       }
       $tags .= '</span> ';
    }

    my $creation = ($exp->creator_id ? $exp->creator->display_name  . ' ' : '') . ($exp->date ne '0000-00-00 00:00:00' ? $exp->date : '');

    my $link;
    $link = "<a target='_blank' href='" . ( $exp->link =~ /^http/ ? $exp->link : 'http://' . $exp->link ) . "'>" . $exp->link . '</a>' if $link;
    my $fields = [
        { title => "ID", value => $exp->id },
        { title => "Name", value => $exp->name},
        { title => "Description", value => $exp->description},
        { title => "Data Type", value => ucfirst($exp->data_type_desc) },
        { title => "Genome", value => $exp->genome->info_html },
        { title => "Source", value => $exp->source->info_html },
        { title => "Link", value => $link },
        { title => "Version", value => $exp->version },
        { title => "Tags", value => $tags || '' },
        { title => "Notebooks", value => $exp->notebooks_desc($EMBED) },
        { title => "Restricted", value => $exp->restricted ? "Yes" : "No"},
        { title => "Creation", value => $creation}
    ];

    my $owner = $exp->owner;
    push @$fields, { title => "Owner", value => $owner->display_name } if $owner;
    
    my $users = ( $exp->restricted ? join(', ', sort map { $_->display_name } $USER->users_with_access($exp)) : 'Everyone' );
    push @$fields, { title => "Users with access", value => $users } if $users;
    
    my $groups = ($exp->restricted ? join(', ', sort map { $_->name } $USER->groups_with_access($exp)) : undef);
    push @$fields, { title => "Groups with access", value => $groups } if $groups;
    
    push @$fields, { title => "Note", value => "This experiment has been deleted" } if $exp->deleted;

    return {
        fields => $fields,
        genome_view_url  => "GenomeView.pl?embed=$EMBED&gid=$gid&tracks=experiment$eid",
        editable => $allow_edit,
        restricted => $exp->restricted
    };
}

sub get_experiment_info {
    my $data = _get_experiment_info(@_);

    my $template_file = catfile($P->{TMPLDIR}, "widgets", "experiment_info_table.tmpl");
    my $info_table = HTML::Template->new(filename => $template_file);

    $info_table->param(fields => $data->{fields});
    $info_table->param(genome_view_url => $data->{genome_view_url});
    $info_table->param(editable => $data->{editable});
    $info_table->param(restricted => $data->{restricted});

    return $info_table->output;
}

sub search_annotation_types {
    my %opts        = @_;
    my $type_group  = $opts{type_group};
    my $search_term = $opts{search_term};

    #print STDERR "search_annotation_types: $search_term $type_group\n";
    return '' unless $search_term;

    $search_term = '%' . $search_term . '%';

    my $group;
    if ($type_group) {
        $group =
          $DB->resultset('AnnotationTypeGroup')->find( { name => $type_group } );
    }

    my @types;
    if ($group) {
        #print STDERR "type_group=$type_group " . $group->id . "\n";
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

sub get_progress_log {
    my %opts = @_;
    my $workflow_id = $opts{workflow_id};
    return unless $workflow_id;

    my (undef, $results_path) = get_workflow_paths($USER->name, $workflow_id);
    return unless (-r $results_path);

    my $result_file = catfile($results_path, '1');
    return unless (-r $result_file);

    my $result = CoGe::Accessory::TDS::read($result_file);
    return unless $result;

    return encode_json(
        {
            experiment_id => $result->{experiment_id},
        }
    );
}

sub get_load_log {
    my %opts = @_;
    my $eid = $opts{eid};
    my $getEverything = $opts{get_everything};
    return unless $eid;
    
    my $log = get_log( item_id => $eid, item_type => 'experiment', getEverything => $getEverything, html => 1 );
    
    return $log;
}

sub send_error_report {
    my %opts = @_;
    my $load_id = $opts{load_id};
    my $job_id = $opts{job_id};

    # Get the staging directory
    my ($staging_dir, $result_dir) = get_workflow_paths($USER->name, $job_id);

    my $url = $P->{SERVER} . "$PAGE_TITLE.pl?";
    $url .= "job_id=$job_id;" if $job_id;
    $url .= "load_id=$load_id" if $load_id;

    my $email = $P->{SUPPORT_EMAIL};

    my $body =
        "SNP finder failed\n\n"
        . 'For user: '
        . $USER->name . ' id='
        . $USER->id . ' '
        . $USER->date . "\n\n"
        . "staging_directory: $staging_dir\n\n"
        . "tiny link: $url\n\n";
    $body .= get_progress_log();

    CoGe::Accessory::Web::send_email(
        from    => $email,
        to      => $email,
        subject => "Load error notification from $PAGE_TITLE",
        body    => $body
    );
}
