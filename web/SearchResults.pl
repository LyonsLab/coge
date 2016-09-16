#! /usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CoGe::Accessory::Web;

use vars qw($CONF $USER $DB %FUNCTION $FORM $SEARCH_TERM $LINK $PAGE_TITLE $PAGE_NAME);

$PAGE_TITLE = 'Search';
$PAGE_NAME  = 'SearchResults.pl';

$FORM = new CGI;
$SEARCH_TERM = $FORM->param('s');

( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init( cgi => $FORM, page_title => $PAGE_TITLE );

%FUNCTION = (
	user_is_admin => \&user_is_admin,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
	my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( USER        => $USER->display_name || '',
	                  PAGE_TITLE  => "Search Results",
	                  TITLE       => "Search Results",
	                  PAGE_LINK   => $LINK,
	                  SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
	                  HOME        => $CONF->{SERVER},
                      HELP        => '',
                      WIKI_URL    => $CONF->{WIKI_URL} || '',
                      ADMIN_ONLY  => $USER->is_admin,
                      CAS_URL     => $CONF->{CAS_URL} || '',
                      COOKIE_NAME => $CONF->{COOKIE_NAME} || '',
                      SEARCH_TERM => $SEARCH_TERM );
	$template->param( LOGON       => 1 ) unless $USER->user_name eq "public";
	$template->param( BODY        => gen_body() );
	
	return $template->output;
}

sub gen_body {
	my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'SearchResults.tmpl' );
	$template->param( API_BASE_URL => $CONF->{SERVER} . 'api/v1/',
	                  USER_NAME   => $USER->user_name,
	                  SEARCH_TERM => $SEARCH_TERM );
	
	return $template->output;
}

sub user_is_admin {
	return $USER->is_admin;
}