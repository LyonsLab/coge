#! /usr/bin/perl -w

use strict;
use CGI;
use HTML::Template;
use CoGeX;
use CoGe::Accessory::Web;

use vars
  qw( $P $PAGE_TITLE $USER $coge %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS
     $EMBED );

$PAGE_TITLE = 'GenomeView';

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

%FUNCTION = ();
CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {    
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param(
            #HELP => "/wiki/index.php?title=$PAGE_TITLE",
	        HELP => $P->{SERVER},
            USER => ( $USER->user_name eq "public" ? '' : $USER->display_name ),
            PAGE_TITLE => 'Genome Viewer',
	        TITLE      => 'Genome Viewer',
            LOGO_PNG   => "CoGe.svg",
            ADJUST_BOX => 1,
            ADMIN_ONLY => $USER->is_admin
        );
        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    }
    
    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );

	my $gid = $FORM->param('gid');
	my $genome = $coge->resultset('Genome')->find($gid);
	return 'Genome not found' unless $genome;
    return 'Access denied' unless ( $USER->has_access_to_genome($genome) );

	$template->param( GENOME_ID => $gid );
	$template->param( GENOME_INFO => $genome->info ) unless $EMBED;
	$template->param( HEIGHT => ($EMBED ? '99%' : '80%') );
	$template->param( WIDTH => ($EMBED ? '99%' : '100%') );

    return $template->output;
}
