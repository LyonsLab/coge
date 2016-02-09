#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadGenome.pl, replicate changes when applicable.

use strict;

use CGI;
use HTML::Template;
use JSON::XS;
use URI::Escape::JavaScript qw(unescape);
#use File::Path;
#use File::Copy;
#use File::Basename;
use File::Spec::Functions qw(catdir catfile);
#use File::Listing qw(parse_dir);
#use File::Slurp;
use LWP::Simple;
use URI;
use Sort::Versions;
use Data::Dumper;

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::IRODS;
use CoGe::Accessory::Jex;
use CoGe::Accessory::TDS;
use CoGe::Accessory::Utils;
use CoGe::Core::Genome qw(genomecmp);
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path get_irods_file);

no warnings 'redefine';

use vars qw(
  $CONF $PAGE_TITLE $TEMPDIR $USER $DB $FORM $LINK $EMBED
  %FUNCTION $LOAD_ID $WORKFLOW_ID $JEX
);

$PAGE_TITLE = 'SynMap3D';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$JEX = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );

# Get workflow_id and load_id for previous load if specified.  Otherwise
# generate a new load_id for data upload.
$WORKFLOW_ID = $FORM->Vars->{'wid'} || $FORM->Vars->{'job_id'}; # wid is new name, job_id is legacy name
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
# AKB TODO $TEMPDIR = get_upload_path($USER->name, $LOAD_ID);

$EMBED = $FORM->param('embed');

%FUNCTION = (
    dotplot_dots        => \&dotplot_dots,
    send_error_report   => \&send_error_report
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub generate_html {
   
    my $template;
    
    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param( 
            PAGE_TITLE   => $PAGE_TITLE,
            TITLE        => "SynMap 3D",
            PAGE_LINK    => $LINK,
            HOME         => $CONF->{SERVER},
            HELP         => 'SynMap3d',
            WIKI_URL     => $CONF->{WIKI_URL} || '',
			ADJUST_BOX   => 1,
            ADMIN_ONLY   => $USER->is_admin,
            USER         => $USER->display_name || '',
            CAS_URL      => $CONF->{CAS_URL} || ''
        );
        $template->param( LOGON      => 1 ) unless $USER->is_public;
    }

    $template->param( BODY => generate_body() );
    return $template->output;
}

sub generate_body {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( 
            PAGE_NAME => "$PAGE_TITLE.pl",
            API_BASE_URL => 'api/v1/'
    );
    
    # Force login 
    if ( $USER->is_public ) {
        $template->param( LOGIN => 1 );
        return $template->output;
    }
    
    # Set genome IDs if specified
    my $x_gid = $FORM->param('x_gid');
    my $y_gid = $FORM->param('y_gid');
    my $z_gid = $FORM->param('z_gid');
    if ($x_gid) {
        my $genome = $DB->resultset('Genome')->find($x_gid);
        if ($genome && $USER->has_access_to_genome($genome)) { # check permission
            $template->param(
                X_GENOME_NAME => $genome->info,
                X_GENOME_ID   => $genome->id
            );
        }
    }
    if ($y_gid) {
        my $genome = $DB->resultset('Genome')->find($y_gid);
        if ($genome && $USER->has_access_to_genome($genome)) { # check permission
            $template->param(
                Y_GENOME_NAME => $genome->info,
                Y_GENOME_ID   => $genome->id
            );
        }
    }
    if ($z_gid) {
        my $genome = $DB->resultset('Genome')->find($z_gid);
        if ($genome && $USER->has_access_to_genome($genome)) { # check permission
            $template->param(
                Z_GENOME_NAME => $genome->info,
                Z_GENOME_ID   => $genome->id
            );
        }
    }

    # Set options if specified
    my $hide_nosynt = $FORM->param('hide');
    my $min_len = $FORM->param('min_len');
    my $sortby = $FORM->param('sortby');
    my $vr = $FORM->param('vr');
    if ($hide_nosynt) {
	$template->param( HIDE_NOSYNT => $hide_nosynt );
    }
    if ($min_len) {}
    if ($sortby) {}
    if ($vr) {}


    $template->param(
        MAIN          => 1,
        PAGE_TITLE    => $PAGE_TITLE,
        EMBED         => $EMBED,
    	LOAD_ID       => $LOAD_ID,
    	WORKFLOW_ID   => $WORKFLOW_ID,
        API_BASE_URL  => 'api/v1/', #TODO move into config file or module
        HELP_URL      => 'https://genomevolution.org/wiki/index.php/SynMap3d',
        SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
        DEFAULT_TAB              => 0,
        USER                     => $USER->user_name
    );
    $template->param( SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
                      SPLASH_CONTENTS    => 'This page allows you to compare synteny between three genomes in a 3d environment. A VR mode is available, for those who so desire.' );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

sub get_debug_log {
    my %opts         = @_;
    my $workflow_id = $opts{workflow_id};
    return unless $workflow_id;
    #TODO authenticate user access to workflow

    my (undef, $results_path) = get_workflow_paths($USER->name, $workflow_id);
    return unless (-r $results_path);

    my $result_file = catfile($results_path, 'debug.log');
    return unless (-r $result_file);

    my $result = read_file($result_file);
    return $result;
}

sub dotplot_dots {
    my %opts = @_;
    my $genome_idX = $opts{genome_idX};
    my $genome_idY = $opts{genome_idY};
    my $genome_idZ = $opts{genome_idZ};
    my $ksfile_xy = $opts{ksfile_xy};
    my $ksfile_xz = $opts{ksfile_xz};
    my $ksfile_yz = $opts{ksfile_yz};
    my $option_name = $opts{option_name};
    my $hide = $opts{hide};
    my $min_len = $opts{min_len};

    my $DIAGSDIR = $CONF->{DIAGSDIR};
    my $SYN3DIR = $CONF->{SYN3DIR};
    my $SCRIPTDIR = $CONF->{SCRIPTDIR};
    my $DIAGSURL = '/asherkhb/coge/data/diags';
    my $SYN3DURL = '/asherkhb/coge/data/syn3d';

    my $dotlog = $option_name . 'log.json';
    my $merge_log = $genome_idX . '_' . $genome_idY . '_' . $genome_idZ . '_' . $option_name . 'log.json';

    my ( $dir1, $dir2 ) = sort ( $genome_idX, $genome_idY );
    my ( $dir3, $dir4 ) = sort ( $genome_idX, $genome_idZ );
    my ( $dir5, $dir6 ) = sort ( $genome_idY, $genome_idZ );

    my $workflow = $JEX->create_workflow( name => "Finding Syntenic Points", init => 1 );

    # Build/add dotplot_dots XY job.
    my $cmd_xy = catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile_xy . ' ' . $option_name;
    my $dots_xy = $dir1 . '_' . $dir2 . '_' . $option_name . 'synteny.json';
    my $dots_xy_log = catfile($DIAGSDIR, $dir1, $dir2, $dir1 . '_' . $dir2 . '_' . $dotlog);
    my $dots_xy_path = catfile($DIAGSDIR, $dir1, $dir2, $dots_xy);
    my $dots_xy_url = catfile($DIAGSURL, $dir1, $dir2, $dots_xy);
    my $outputs_xy = [$dots_xy_log, $dots_xy_path];
    $workflow->add_job({
        cmd => $cmd_xy,
        outputs => $outputs_xy,
        description => "running XY dotplot_dots...",
    });

    # Build/add dotplot_dots XZ job.
    my $cmd_xz = catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile_xz . ' ' . $option_name;
    my $dots_xz = $dir3 . '_' . $dir4  . '_' . $option_name . 'synteny.json';
    my $dots_xz_log = catfile($DIAGSDIR, $dir3, $dir4, $dir3 . '_' . $dir4 . '_' . $dotlog);
    my $dots_xz_path = catfile($DIAGSDIR, $dir3, $dir4, $dots_xz);
    my $dots_xz_url = catfile($DIAGSURL, $dir3, $dir4, $dots_xz);
    my $outputs_xz = [$dots_xz_log, $dots_xz_path];
    $workflow->add_job({
        cmd => $cmd_xz,
        outputs => $outputs_xz,
        description => "running XZ dotplot_dots...",
    });

    # Build/add dotplot_dots YZ job.
    my $cmd_yz = catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile_yz . ' ' . $option_name;
    my $dots_yz = $dir5 . '_' . $dir6 . '_' . $option_name . 'synteny.json';
    my $dots_yz_log = catfile($DIAGSDIR, $dir5, $dir6, $dir5 . '_' . $dir6 . '_' . $dotlog);
    my $dots_yz_path = catfile($DIAGSDIR, $dir5, $dir6, $dots_yz);
    my $dots_yz_url = catfile($DIAGSURL, $dir5, $dir6, $dots_yz);
    my $outputs_yz = [$dots_yz_log, $dots_yz_path];
    $workflow->add_job({
        cmd => $cmd_yz,
        outputs => $outputs_yz,
        description => "running YZ dotplot_dots...",
    });

    # Build/add three_dots_merge job.
    # my $merge_ins = ' -i1 ' . @$outputs_xy[1] . ' -i2 ' . @$outputs_xz[1] . ' -i3 ' . @$outputs_yz[1];
    my $merge_ins = ' -i1 ' . $dots_xy_path . ' -i2 ' . $dots_xz_path . ' -i3 ' . $dots_yz_path;
    my $merge_ots = ' -o ' . $SYN3DIR . ' -n ' . '"' . $option_name . '"';
    my $merge_gids = ' -xid ' . $genome_idX . ' -yid ' . $genome_idY . ' -zid ' . $genome_idZ;
    my $merge_opts = '';
    if ($hide eq 'true') { $merge_opts .= ' -P' }
    if ($min_len > 0) { $merge_opts .= ' -M ' . $min_len; }

    my $merge_cmd = catfile($SCRIPTDIR, 'three_dots_merge.py') . $merge_ins . $merge_ots . $merge_gids . $merge_opts;
    my $merge_in = [$dots_xy_path, $dots_xz_path, $dots_yz_path];
    my $dots_xyz = $genome_idX . '_' . $genome_idY . '_' . $genome_idZ . '_' . $option_name . 'dots.json';
    my $dots_xyz_path = catfile($SYN3DIR, $dots_xyz);
    my $dots_xyz_url = catfile($SYN3DURL, $dots_xyz);
    my $hist_xyz = $genome_idX . '_' . $genome_idY . '_' . $genome_idZ . '_' . $option_name . 'histogram.json';
    my $hist_xyz_path = catfile($SYN3DIR, $hist_xyz);
    my $hist_xyz_url = catfile($SYN3DURL, $hist_xyz);
    my $merge_out = [$dots_xyz_path, catfile($SYN3DIR, $merge_log), $hist_xyz_path];
    $workflow->add_job({
        cmd => $merge_cmd,
        inputs => $merge_in,
        outputs => $merge_out,
        description => "merging XYZ dots..."
    });

    my $response = $JEX->submit_workflow($workflow);
    return encode_json(
        {
            id => $response->{id},
            status  => $response->{status},
            success => $JEX->is_successful($response)
            ? JSON::true
            : JSON::false,
            xy_json => $dots_xy_url,
            xz_json => $dots_xz_url,
            yz_json => $dots_yz_url,
            merge_json => $dots_xyz_url,
            histo_json => $hist_xyz_url
        }
    );
}

sub send_error_report {
    my %opts = @_;
    my $load_id = $opts{load_id};
    my $job_id = $opts{job_id};
    unless ($load_id and $job_id) {
        print STDERR "LoadExperiment::send_error_report: missing required params\n";
        return;
    }

    # Get the staging directory
    my ($staging_dir, $result_dir) = get_workflow_paths($USER->name, $job_id);

    my $url = $CONF->{SERVER} . "$PAGE_TITLE.pl?";
    $url .= "job_id=$job_id;" if $job_id;
    $url .= "load_id=$load_id";

    my $email = $CONF->{SUPPORT_EMAIL};

    my $body =
        "Load failed\n\n"
        . 'For user: '
        . $USER->name . ' id='
        . $USER->id . ' '
        . $USER->date . "\n\n"
        . "staging_directory: $staging_dir\n\n"
        . "result_directory: $result_dir\n\n"
        . "tiny link: $url\n\n";

    my $log = get_debug_log(workflow_id => $job_id);
    $body .= $log if $log;

    CoGe::Accessory::Web::send_email(
        from    => $email,
        to      => $email,
        subject => "Load error notification from $PAGE_TITLE",
        body    => $body
    );
}
