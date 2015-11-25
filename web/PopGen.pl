#! /usr/bin/perl -w
use strict;
use CoGe::Accessory::Web;
use CoGe::Algos::PopGen::Results qw( export get_data );
use CGI;
use Data::Dumper;
use File::Spec::Functions qw( catfile );
use HTML::Template;
use JSON::XS;
use POSIX;
no warnings 'redefine';

use vars qw($CONF $USER $FORM $DB $PAGE_TITLE $LINK);

$PAGE_TITLE = 'PopGen';

$| = 1;    # turn off buffering

$FORM = new CGI;

my $export = $FORM->Vars->{'export'};
if ($export) {
	export catfile($ENV{COGE_HOME}, 'test.txt'), $FORM->Vars->{'type'}, $export;
	return;
}

( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

print $FORM->header, gen_html();

sub gen_html {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => $PAGE_TITLE,
		              TITLE      => $PAGE_TITLE,
    				  PAGE_LINK  => $LINK,
				      HOME       => $CONF->{SERVER},
                      HELP       => $PAGE_TITLE,
                      WIKI_URL   => $CONF->{WIKI_URL} || '',
                      ADJUST_BOX => 1,
    				  ADMIN_ONLY => $USER->is_admin,
                      CAS_URL    => $CONF->{CAS_URL} || '',
                      USER       => $USER->display_name || ''
    );
    $template->param( LOGON => 1 ) unless ($USER->user_name eq "public");
    $template->param( BODY => gen_body() );

    return $template->output;
}

sub gen_body {
    if (!$FORM->Vars->{'id'}) {
        return 'id url parameter not set';
    }
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    return $template->output;
}
