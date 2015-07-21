#! /usr/bin/perl -w

use warnings;
use strict;

use CGI;
use HTML::Template;
use JSON::XS;
use URI::Escape::JavaScript qw(unescape);
use File::Path;
use File::Copy;
use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use File::Listing qw(parse_dir);
use Sort::Versions;
use Data::Dumper;

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::TDS;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path);

use vars qw(
  $CONF $PAGE_TITLE $TEMPDIR $USER $DB $FORM $LINK $EMBED
  %FUNCTION $LOAD_ID $WORKFLOW_ID
);

$PAGE_TITLE = 'Project';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

# Get workflow_id and load_id for previous load if specified.  Otherwise
# generate a new load_id for data upload.
$WORKFLOW_ID = $FORM->Vars->{'wid'} || $FORM->Vars->{'job_id'}; # wid is new name, job_id is legacy name
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = get_upload_path($USER->name, $LOAD_ID);

$EMBED = $FORM->param('embed');

%FUNCTION = (
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
        $template->param( PAGE_TITLE => $PAGE_TITLE,
		                  TITLE      => "Project",
        	              PAGE_LINK  => $LINK,
			              HOME       => $CONF->{SERVER},
                          HELP       => 'Project',
                          WIKI_URL   => $CONF->{WIKI_URL} || '',
			              ADJUST_BOX => 1,
                          ADMIN_ONLY => $USER->is_admin,
                          USER       => $USER->display_name || '',
                          CAS_URL    => $CONF->{CAS_URL} || ''
        );
        $template->param( LOGON      => 1 ) unless $USER->is_public;
    }

    $template->param( BODY => generate_body() );
    return $template->output;
}

sub generate_body {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
    
    # Force login
    if ( $USER->is_public ) {
        $template->param( LOGIN => 1 );
        return $template->output;
    }

    # Set genome ID if specified
    my $gid = $FORM->param('gid');
    if ($gid) {
        my $genome = $DB->resultset('Genome')->find($gid);
        if ($genome && $USER->has_access_to_genome($genome)) { # check permission
            $template->param(
                GENOME_NAME => $genome->info,
                GENOME_ID   => $genome->id
            );
        }
    }
    
    $template->param(
        MAIN          => 1,
        PAGE_TITLE    => $PAGE_TITLE,
        EMBED         => $EMBED,
    	LOAD_ID       => $LOAD_ID,
    	WORKFLOW_ID   => $WORKFLOW_ID,
        API_BASE_URL  => 'api/v1/', #TODO move into config file or module
        HELP_URL      => 'https://genomevolution.org/wiki/index.php/Project',
        SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
        USER          => $USER->user_name
    );
#    $template->param( SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
#                      SPLASH_CONTENTS    => 'This page allows you to load quantitative, polymorphism, or alignment data onto a genome from a variety of file formats.' );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}
