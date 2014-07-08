#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use CoGe::Accessory::IRODS;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage qw(get_workflow_paths);
use HTML::Template;
use JSON::XS;
use Spreadsheet::WriteExcel;
use File::Basename;
use File::Path;
use File::Slurp;
use File::Spec::Functions qw(catdir catfile);
use Sort::Versions;
use CoGe::Pipelines::FindSNPs qw( run );
use Data::Dumper;

use vars qw(
    $P $PAGE_TITLE $USER $LINK $coge $FORM $EMBED %FUNCTION $ERROR
    $JOB_ID $LOAD_ID $TEMPDIR $CONFIGFILE
);

$PAGE_TITLE = "ExperimentView";
$ERROR = encode_json( { error => 1 } );
$CONFIGFILE = $ENV{COGE_HOME} . '/coge.conf';

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE,
);

$JOB_ID  = $FORM->Vars->{'job_id'};
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = $P->{SECTEMPDIR} . $PAGE_TITLE . '/' . $USER->name . '/' . $LOAD_ID . '/';

%FUNCTION = (
    get_experiment_info        => \&get_experiment_info,
    edit_experiment_info       => \&edit_experiment_info,
    update_experiment_info     => \&update_experiment_info,
    get_sources                => \&get_sources,
    make_experiment_public     => \&make_experiment_public,
    make_experiment_private    => \&make_experiment_private,
    add_type_to_experiment     => \&add_type_to_experiment,
    get_experiment_types       => \&get_experiment_types,
    get_type_description       => \&get_type_description,
    remove_experiment_type     => \&remove_experiment_type,
    get_annotations            => \&get_annotations,
    add_annotation             => \&add_annotation,
    update_annotation          => \&update_annotation,
    remove_annotation          => \&remove_annotation,
    get_annotation             => \&get_annotation,
    search_annotation_types    => \&search_annotation_types,
    get_annotation_type_groups => \&get_annotation_type_groups,
    check_login                => \&check_login,
    export_experiment_irods    => \&export_experiment_irods,
    get_file_urls              => \&get_file_urls,
    find_snps                  => \&find_snps,
    get_progress_log           => \&get_progress_log,
    send_error_report          => \&send_error_report,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub edit_experiment_info {
    my %opts = @_;
    my $eid  = $opts{eid};
    return 0 unless $eid;

    my $exp = $coge->resultset('Experiment')->find($eid);
    my $desc = ( $exp->description ? $exp->description : '' );

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        EDIT_EXPERIMENT_INFO => 1,
        EID                  => $eid,
        NAME                 => $exp->name,
        DESC                 => $desc,
        SOURCE               => $exp->source->name,
        SOURCE_ID            => $exp->source->id,
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
    my $version   = $opts{version};

    my $exp = $coge->resultset('Experiment')->find($eid);
    $exp->name($name);
    $exp->description($desc) if $desc;
    $exp->version($version);
    $exp->data_source_id($source_id);
    $exp->update;

    return 1;
}

sub get_sources {
    #my %opts = @_;

    my @sources;
    foreach ( $coge->resultset('DataSource')->all() ) {
        push @sources, { 'label' => $_->name, 'value' => $_->id };
    }

    return encode_json( \@sources );
}

sub make_experiment_public {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No EID specified" unless $eid;
    #return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $exp = $coge->resultset('Experiment')->find($eid);
    $exp->restricted(0);
    $exp->update;

    return 1;
}

sub make_experiment_private {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No EID specified" unless $eid;
    #return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $exp = $coge->resultset('Experiment')->find($eid);
    $exp->restricted(1);
    $exp->update;

    return 1;
}

sub genomecmp {    # FIXME mdb 8/24/12 - redundant declaration in ListView.pl
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    versioncmp( $b->version, $a->version )
      || $a->type->id <=> $b->type->id
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

sub add_type_to_experiment {
    my %opts = @_;
    my $eid  = $opts{eid};
    return 0 unless $eid;
    my $name = $opts{name};
    return 0 unless $name;
    my $description = $opts{description};

    my $type = $coge->resultset('ExperimentType')->find( { name => $name, description => $description } );

    if ($type) {
        # If type exists, check if already assigned to this experiment
        foreach ( $type->experiment_type_connectors ) {
            return 1 if ( $_->experiment_id == $eid );
        }
    }
    else {
        # Create type if it doesn't already exist
        $type =
          $coge->resultset('ExperimentType')
          ->create( { name => $name, description => $description } );
    }

    # Create connection
    $coge->resultset('ExperimentTypeConnector')->create(
        {
            experiment_id      => $eid,
            experiment_type_id => $type->id
        }
    );

    return 1;
}

sub get_experiment_types {
    #my %opts = @_;

    my %unique;

    my $rs = $coge->resultset('ExperimentType');
    while ( my $et = $rs->next ) {
        $unique{ $et->name }++;
    }

    return encode_json( [ sort keys %unique ] );
}

sub get_type_description {
    my %opts = @_;
    my $name = $opts{name};
    return unless ($name);

    my $type = $coge->resultset('ExperimentType')->find( { name => $name } );
    return $type->description if ($type);
}

sub linkify {
    my ( $link, $desc ) = @_;
    return "<span class='link' onclick=\"window.open('$link')\">$desc</span>";
}

sub remove_experiment_type {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "No experiment ID specified" unless $eid;
    my $etid = $opts{etid};
    return "No experiment type ID specified" unless $etid;
    #return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );

    my $etc =
      $coge->resultset('ExperimentTypeConnector')
      ->find( { experiment_id => $eid, experiment_type_id => $etid } );
    $etc->delete();

    return 1;
}

sub get_annotations {
    my %opts = @_;
    my $eid  = $opts{eid};
    return "Must have valid experiment id\n" unless ($eid);

    my $exp = $coge->resultset('Experiment')->find($eid);
    return "Access denied\n" unless $USER->has_access_to_experiment($exp);

    my $user_can_edit = ( $USER->is_admin || $USER->is_owner_editor( experiment => $eid ) );

    # Categorize annotations based on type group and type
    my %groups;
    my $num_annot = 0;
    foreach my $a ( $exp->annotations ) {
        my $group = ( $a->type->group ? $a->type->group->name : undef);
        my $type = $a->type->name;
        push @{ $groups{$group}{$type} }, $a if (defined $group and defined $type);
        $num_annot++;
    }

    # Build annotation table
    my $html;
    if ($num_annot) {
        $html .= '<table id="experiment_annotation_table" class="ui-widget-content ui-corner-all small" style="max-width:800px;overflow:hidden;word-wrap:break-word;border-spacing:0;"><thead style="display:none"></thead><tbody>';
        foreach my $group ( sort keys %groups ) { # groups
            my $first_group = 1;
            foreach my $type ( sort keys %{ $groups{$group} } ) { # types
                my $first_type = 1;
                foreach my $a ( sort { $a->id <=> $b->id } @{ $groups{$group}{$type} } ) { # annotations
                    my $header = ($group and $first_group-- > 0 ? "<b>$group</b>: " : '') . ($first_type-- > 0 ? "$type:" : '');
                    $html .= "<tr style='vertical-align:top;'>";
                    $html .= "<th align='right' class='title5' style='padding-right:10px;white-space:nowrap;font-weight:normal;background-color:white;'>$header</th>";
                    #$html .= '<td>';
                    my $image_link = ( $a->image ? 'image.pl?id=' . $a->image->id : '' );
                    my $image_info = (
                        $a->image
                        ? "<a href='$image_link' target='_blank' title='click for full-size image'><img height='40' width='40' src='$image_link' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='float:left;padding:1px;border:1px solid lightgray;margin-right:5px;'></a>"
                        : ''
                    );
                    #$html .= $image_info if $image_info;
                    #$html .= "</td>";
                    $html .= "<td class='data5'>" . $image_info . $a->info . '</td>';
                    $html .= '<td style="padding-left:5px;">';
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
        $html .= '<table class="ui-widget-content ui-corner-all small padded note"><tr><td>There a no additional metadata items for this experiment.</tr></td></table>';
    }

    if ($user_can_edit) {
        $html .= qq{<span onClick="add_annotation_dialog();" style="font-size: .75em" class='ui-button ui-button-icon-left ui-corner-all'><span class="ui-icon ui-icon-plus"></span>Add</span>};
    }

    return $html;
}

sub get_annotation {
    my %opts = @_;
    my $aid  = $opts{aid};
    return unless $aid;

    #TODO check user access here

    my $ea = $coge->resultset('ExperimentAnnotation')->find($aid);
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
    my $eid  = $opts{parent_id};
    return 0 unless $eid;
    my $type_group = $opts{type_group};
    my $type       = $opts{type};
    return 0 unless $type;
    my $annotation     = $opts{annotation};
    my $link           = $opts{link};
    my $image_filename = $opts{edit_annotation_image};
    my $fh             = $FORM->upload('edit_annotation_image');

    #print STDERR "add_annotation: $eid $type $annotation $link\n";

    # Create the type and type group if not already present
    my $group_rs;
    if ($type_group) {
        $group_rs =
          $coge->resultset('AnnotationTypeGroup')
          ->find_or_create( { name => $type_group } );
    }

    my $type_rs = $coge->resultset('AnnotationType')->find_or_create(
        {
            name                     => $type,
            annotation_type_group_id => ( $group_rs ? $group_rs->id : undef )
        }
    );

    # Create the image
    my $image;
    if ($fh) {
        read( $fh, my $contents, -s $fh );
        $image = $coge->resultset('Image')->create(
            {
                filename => $image_filename,
                image    => $contents
            }
        );
        return 0 unless $image;
    }

    # Create the annotation
    my $annot = $coge->resultset('ExperimentAnnotation')->create(
        {
            experiment_id      => $eid,
            annotation         => $annotation,
            link               => $link,
            annotation_type_id => $type_rs->id,
            image_id           => ( $image ? $image->id : undef )
        }
    );
    return 0 unless $annot;

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

    my $ea = $coge->resultset('ExperimentAnnotation')->find($aid);
    return unless $ea;

    # Create the type and type group if not already present
    my $group_rs;
    if ($type_group) {
        $group_rs =
          $coge->resultset('AnnotationTypeGroup')
          ->find_or_create( { name => $type_group } );
    }
    my $type_rs = $coge->resultset('AnnotationType')->find_or_create(
        {
            name                     => $type,
            annotation_type_group_id => ( $group_rs ? $group_rs->id : undef )
        }
    );

    # Create the image
    #TODO if image was changed delete previous image
    my $image;
    if ($fh) {
        read( $fh, my $contents, -s $fh );
        $image = $coge->resultset('Image')->create(
            {
                filename => $image_filename,
                image    => $contents
            }
        );
        return 0 unless $image;
    }

    $ea->annotation($annotation);
    $ea->link($link);
    $ea->annotation_type_id( $type_rs->id );
    $ea->image_id( $image->id ) if ($image);
    $ea->update;

    return;
}

#XXX: Move to a module
sub check_login {
    #print STDERR $USER->user_name . ' ' . int($USER->is_public) . "\n";
    return ($USER && !$USER->is_public);
}

#XXX: Move to a module
sub get_irods_path {
    my $username = $USER->user_name;
    my $dest = $P->{IRODSDIR};
    $dest =~ s/\<USER\>/$username/;
    return $dest;
}

sub export_experiment_irods {
    my %opts = @_;
    my $eid = $opts{eid};

    my $experiment = $coge->resultset('Experiment')->find($eid);
    return $ERROR unless $USER->has_access_to_experiment($experiment);

    my ($statusCode, $file) = generate_export($eid);

    unless($statusCode) {
        my $genome = $experiment->genome;
        my @types = $experiment->types;
        my @notebooks = $experiment->notebooks;
        my $dir = get_irods_path();
        my $dest = File::Spec->catdir($dir, basename($file));
        my $restricted = ($experiment->restricted) ? "yes" : "no";

        my %meta = (
            'Imported From'            => "CoGe: " . $P->{SERVER},
            'CoGe ExperimentView Link' => $P->{SERVER} . "ExperimentView.pl?eid=$eid",
            'CoGe GenomeInfo Link'     => $P->{SERVER} . "GenomeInfo.pl?gid=" . $genome->id,
            'Source'                   => $experiment->source->info,
            'Version'                  => $experiment->version,
            'Restricted'               => $restricted
        );

        my $genome_name = $genome->info;
        $genome_name =~ s/&reg;\s*//;

        $meta{'Name'} = $experiment->name if ($experiment->name);
        $meta{'Description'} = $experiment->description if ($experiment->description);
        $meta{'Genome'} = $genome_name;
        $meta{'Source Link'} = $experiment->source->link if $experiment->source->link;
        $meta{'Rows'} = commify($experiment->row_count);

        my $i = 1;
        foreach my $type (@types) {
            my $key = (scalar @types > 1) ? "Experiment Type $i" : "Experiment Type";
            $meta{$key} = $type->name;
            $i++;
        }

        $i = 1;
        foreach my $type (@notebooks) {
            my $key = (scalar @notebooks > 1) ? "Notebook $i" : "Notebook";
            $meta{$key} = $type->name;
            $i++;
        }

        foreach my $a ( $experiment->annotations ) {
            my $group = (
                defined $a->type->group
                ? $a->type->group->name . ',' . $a->type->name
                : $a->type->name
            );

            $meta{$group} = $a->info;
        }

        CoGe::Accessory::IRODS::irods_iput($file, $dest);
        CoGe::Accessory::IRODS::irods_imeta($dest, \%meta);
    }

    return basename($file);
}

sub generate_export {
    my $eid = shift;
    my $filename = "experiment_$eid.tar.gz";

    my $conf = File::Spec->catdir($P->{COGEDIR}, "coge.conf");
    my $script = File::Spec->catdir($P->{SCRIPTDIR}, "export_experiment.pl");
    my $workdir = get_download_path($eid);
    my $resdir = $P->{RESOURCESDIR};

    my $cmd = "$script -eid $eid -config $conf -dir $workdir -output $filename -a 1";

    return (execute($cmd),  File::Spec->catdir(($workdir, $filename)));
}

sub get_download_path {
    my $unique_path = get_unique_id();
    my @paths = ($P->{SECTEMPDIR}, "ExperimentView/downloads", shift, $unique_path);
    return File::Spec->catdir(@paths);
}

sub get_download_url {
    my %args = @_;
    my $id = $args{id};
    my $dir = $args{dir};
    my $filename = basename($args{file});

    my @url = ($P->{SERVER}, "services/JBrowse",
        "service.pl/download/ExperimentView",
        "?eid=$id&dir=$dir&file=$filename");

    return join "/", @url;
}

sub get_file_urls {
    my %opts = @_;
    my $eid = $opts{eid};

    my $experiment = $coge->resultset('Experiment')->find($eid);
    return 0 unless $USER->has_access_to_experiment($experiment);

    my ($statusCode, $file) = generate_export($eid);

    unless($statusCode) {
        my $dir = basename(dirname($file));
        my $url = get_download_url(id => $eid, dir => $dir, file => $file);
        return encode_json({ files => [$url] });
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

    my $ea = $coge->resultset('ExperimentAnnotation')->find( { experiment_annotation_id => $eaid } );
    $ea->delete();

    return 1;
}

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template =
          HTML::Template->new(
            filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template =
          HTML::Template->new(
            filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= ' ' . $USER->last_name
          if ( $USER->first_name && $USER->last_name );
        $template->param(
            PAGE_TITLE => $PAGE_TITLE,
            PAGE_LINK  => $LINK,
            HELP       => '/wiki/index.php?title=' . $PAGE_TITLE,
            USER       => $name,
            LOGO_PNG   => "$PAGE_TITLE-logo.png",
            ADJUST_BOX => 1,
        );
        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    }

    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $eid = $FORM->param('eid');
    return "Need a valid experiment id\n" unless $eid;

    my $exp = $coge->resultset('Experiment')->find($eid);
    return "Access denied" unless $USER->has_access_to_experiment($exp);

    my $gid = $exp->genome_id;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        MAIN            => 1,
        PAGE_NAME       => $PAGE_TITLE . '.pl',
        EID             => $eid,
        DEFAULT_TYPE    => 'note',
        rows            => commify($exp->row_count),
        IRODS_HOME      => get_irods_path(),
        JOB_ID          => $JOB_ID,
        STATUS_URL      => 'jex/status/',
        ALIGNMENT_TYPE  => ($exp->data_type == 3), # FIXME: hardcoded type value
        PUBLIC          => $USER->user_name eq "public" ? 1 : 0
    );
    $template->param( EXPERIMENT_INFO => get_experiment_info( eid => $eid ) || undef );
    $template->param( EXPERIMENT_ANNOTATIONS => get_annotations( eid => $eid ) || undef );

    return $template->output;
}

sub get_experiment_info {
    my %opts  = @_;
    my $eid   = $opts{eid};
    my ($exp) = $coge->resultset('Experiment')->find($eid);
    return "Access denied\n" unless $USER->has_access_to_experiment($exp);

    return "Unable to find an entry for $eid" unless $exp;

    my $allow_edit = $USER->is_admin || $USER->is_owner_editor( experiment => $eid );

    my $gid = $exp->genome->id;

    my $html;
    $html .= $exp->annotation_pretty_print_html( allow_delete => $allow_edit );
    $html .= qq{<a style="font-size: .75em; color: black; float:right;" target="_blank" class='ui-button ui-corner-all ui-button-icon-right' href="GenomeView.pl?gid=$gid&tracks=experiment$eid">View<span class="ui-icon ui-icon-extlink"></span></a>};

    $html .= "<div class='inline'>";

    if ($allow_edit) {
        $html .= qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="edit_experiment_info();">Edit Info</span>};
        $html .= qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="\$('#experiment_type_edit_box').dialog('open');">Add Type</span>};
    }

    if ( $USER->is_admin || $USER->is_owner( experiment => $eid ) ) {
        if ( $exp->restricted ) {
            $html .= qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="make_experiment_public();">Make Public</span>};
        }
        else {
            $html .= qq{<span style="font-size: .75em" class='ui-button ui-corner-all' onClick="make_experiment_private();">Make Private</span>};
        }
    }

    $html .= "</div>";

    return $html;
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
          $coge->resultset('AnnotationTypeGroup')->find( { name => $type_group } );
    }

    my @types;
    if ($group) {
        #print STDERR "type_group=$type_group " . $group->id . "\n";
        @types = $coge->resultset("AnnotationType")->search(
            \[
                'annotation_type_group_id = ? AND (name LIKE ? OR description LIKE ?)',
                [ 'annotation_type_group_id', $group->id ],
                [ 'name',                     $search_term ],
                [ 'description',              $search_term ]
            ]
        );
    }
    else {
        @types = $coge->resultset("AnnotationType")->search(
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
    my %unique;

    my $rs = $coge->resultset('AnnotationTypeGroup');
    while ( my $atg = $rs->next ) {
        $unique{ $atg->name }++;
    }

    return encode_json( [ sort keys %unique ] );
}

sub find_snps {
    my %opts        = @_;
    my $user_name   = $opts{user_name};
    my $eid         = $opts{eid};

    # Check login
    if ( !$user_name || !$USER->is_admin ) {
        $user_name = $USER->user_name;
    }
    if ($user_name eq 'public') {
        return encode_json({ error => 'Not logged in' });
    }

    # Get experiment
    my $experiment = $coge->resultset('Experiment')->find($eid);
    return encode_json({ error => 'Experiment not found' }) unless $experiment;

    # Submit workflow to generate experiment
    my ($workflow_id, $error_msg) = CoGe::Pipelines::FindSNPs::run(
        db => $coge,
        experiment => $experiment,
        user => $USER
    );
    unless ($workflow_id) {
        print STDERR $error_msg, "\n";
        return encode_json({ error => "Workflow submission failed: " . $error_msg });
    }

    # Get tiny link
    my $link = CoGe::Accessory::Web::get_tiny_link(
        url => $P->{SERVER} . "$PAGE_TITLE.pl?job_id=" . $workflow_id
    );

    return encode_json({ job_id => $workflow_id, link => $link });
}

sub get_progress_log {
    my %opts         = @_;
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
#sub get_progress_log {
#    my $logfile = catfile($TEMPDIR, 'staging', 'load_experiment', 'log.txt');
#    open( my $fh, $logfile ) or
#        return encode_json( { status => -1, log => "Error opening log file" } );
#
#    my @lines = ();
#    my ($eid, $nid, $new_load_id);
#    my $status = 0;
#    my $message = '';
#    while (<$fh>) {
#        push @lines, $1 if ( $_ =~ /^log: (.+)/i );
#        if ( $_ =~ /All done/i ) {
#            $status = 1;
#
#            # Generate a new load session ID in case the user chooses to
#            # reuse the form to start another load.
#            $new_load_id = get_unique_id();
#
#            last;
#        }
#        elsif ( $_ =~ /experiment id: (\d+)/i ) {
#            $eid = $1;
#        }
#        elsif ( $_ =~ /log: error: input file is empty/i ) {
#            $status = -2;
#            $message = 'No SNPs were detected in this experiment';
#            last;
#        }
#        elsif ( $_ =~ /log: error/i ) {
#            $status = -1;
#            last;
#        }
#    }
#
#    close($fh);
#
#    return encode_json(
#        {
#            status        => $status,
#            experiment_id => $eid,
#            new_load_id   => $new_load_id,
#            message       => $message
#        }
#    );
#}

sub send_error_report {
    my %opts = @_;
    my $load_id = $opts{load_id};
    my $job_id = $opts{job_id};

    my $staging_dir = catdir($TEMPDIR, 'staging');

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
