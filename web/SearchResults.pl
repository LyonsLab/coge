#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(format_time_diff);
use CoGe::Accessory::Jex;

$|=1;

use vars
  qw($P $PAGE_NAME $USER $BASEFILE $coge $cogeweb %FUNCTION $FORM $MAX_SEARCH_RESULTS %ITEM_TYPE $JEX);

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init( cgi => $FORM );

$JEX =
  CoGe::Accessory::Jex->new( host => $P->{JOBSERVER}, port => $P->{JOBPORT} );

%FUNCTION = (
	user_is_admin					=> \&user_is_admin,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
	#print STDERR "HTML\n";
	my $html;
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( USER       => $USER->display_name || '',
	                  PAGE_TITLE => qq{Search Results},
	                  TITLE      => "Search Results",
	                  HOME       => $P->{SERVER},
                      HELP       => '',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      CAS_URL    => $P->{CAS_URL} || '',
                      ADMIN_ONLY => $USER->is_admin );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( BODY       => gen_body() );
	$html .= $template->output;
}

sub gen_body {
	
	my $term = $FORM->param('s');
	say STDERR "$term\n";

	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'SearchResults.tmpl' );
	$template->param( 	MAIN 			=> 1,
						API_BASE_URL  	=> 'api/v1/',  
						USER			=> $USER->user_name,
						TERM			=> $term,	
					);
					
	return $template->output;
}

sub user_is_admin {
	return $USER->is_admin;
}