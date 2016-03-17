#! /usr/bin/perl -w
use strict;
use CoGe::Accessory::Web;
use CGI;
use Data::Dumper;
use HTML::Template;
use JSON::XS;
use POSIX;
no warnings 'redefine';

use vars qw($P $USER $FORM $ACCN $FID $db $PAGE_NAME $PAGE_TITLE $LINK);

$PAGE_TITLE = 'My Account';
$PAGE_NAME  = "Account.pl";

$| = 1;    # turn off buffering

$FORM = new CGI;

( $db, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my %FUNCTION = (
);
CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'Account.tmpl' );
    $template->param( PAGE_TITLE => 'My Account',
                      TITLE      => 'View my account info',
                      PAGE_LINK  => $LINK,
                      HEAD       => qq{},
                      HOME       => $P->{SERVER},
                      HELP       => 'My Account',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      ADMIN_ONLY => $USER->is_admin,
                      USER       => $USER->display_name || '',
                      CAS_URL    => $P->{CAS_URL} || ''
    );
    $template->param( LOGON    => 1 ) unless ($USER->user_name eq "public");
    $template->param( email => $USER->email );
    $template->param( image => $USER->image_id );
    $template->param( date => $USER->date );
    #$template->param( description => $USER->description );
    $template->param( username => $USER->user_name );
    return $template->output;
}
