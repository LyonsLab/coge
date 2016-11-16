#! /usr/bin/perl -w

use strict;

use CGI;
use HTML::Template;
use File::Spec::Functions qw(catdir);

use CoGeX;
use CoGe::Accessory::Web;

use vars qw( 
    $CONF $PAGE_TITLE $USER $DB %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS
    $EMBED $LINK
);

$PAGE_TITLE = 'GenomeView';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

%FUNCTION = ();
CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {    
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param(
            USER       => $USER->display_name || '',
            PAGE_TITLE => 'Genome Viewer',
            PAGE_LINK  => $LINK,
            SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
	        TITLE      => 'Genome Viewer',
	        HOME       => $CONF->{SERVER},
            HELP       => 'GenomeView',
            WIKI_URL   => $CONF->{WIKI_URL} || '',
            ADMIN_ONLY => $USER->is_admin,
            CAS_URL    => $CONF->{CAS_URL} || '',
            NO_DOCTYPE => 1,
            COOKIE_NAME => $CONF->{COOKIE_NAME} || ''
        );
        $template->param( LOGON => 1 ) unless ($USER->user_name eq 'public');
    }
    
    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );

    # Get specified genome and verify permissions
	my $gid = $FORM->param('gid');
    if ($gid) {
    	my $genome = $DB->resultset('Genome')->find($gid);
	    return 'Genome not found' unless $genome;
        return 'Access denied' unless ( $USER->has_access_to_genome($genome) );
    	$template->param( GENOME_INFO => $genome->info ) unless $EMBED;
    }
	$template->param( GENOME_ID => $gid,
                      EMBED => $EMBED ? 1 : 0,
	                  HEIGHT => ($EMBED ? '99%' : '80%'),
	                  WIDTH => ($EMBED ? '99%' : '100%'),
	                  API_BASE_URL => $CONF->{SERVER} . 'api/v1/jbrowse', # mdb added base URL, 2/3/15 COGE-289
	                  USER_NAME => $USER->user_name
	);

    return $template->output;
}
