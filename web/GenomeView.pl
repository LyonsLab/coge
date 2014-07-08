#! /usr/bin/perl -w

use strict;
use CGI;
use HTML::Template;
use CoGeX;
use CoGe::Accessory::Web;

use vars
  qw($P $PAGE_TITLE $USER $coge %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS);

$PAGE_TITLE = 'GenomeView';

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

%FUNCTION = ();
CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param(
        HELP => "/wiki/index.php?title=$PAGE_TITLE",
        USER => ( $USER->user_name eq "public" ? '' : $USER->display_name ),
        PAGE_TITLE => 'Genome Viewer',
        LOGO_PNG   => "$PAGE_TITLE-logo.png",
        ADJUST_BOX => 1,
        BODY       => gen_body()
    );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";

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
	$template->param( GENOME_INFO => $genome->info );

    return $template->output;
}
