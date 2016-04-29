#! /usr/bin/perl -w

use strict;

use CGI;
use HTML::Template;
use CoGe::Accessory::Web;

use vars qw(
    $CONF $PAGE_TITLE $USER $LINK $DB $FORM
);

$PAGE_TITLE = 'Error';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    page_title => $PAGE_TITLE,
    cgi => $FORM
);

CoGe::Accessory::Web->dispatch( $FORM, undef, \&gen_html );

sub gen_html {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( USER       => $USER->display_name || '',
                      PAGE_TITLE => 'Error',
				      TITLE      => "Error",
    				  PAGE_LINK  => $LINK,
    				  SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
    				  HOME       => $CONF->{SERVER},
                      HELP       => '',
                      WIKI_URL   => $CONF->{WIKI_URL} || '',
                      ADMIN_ONLY => $USER->is_admin,
                      CAS_URL    => $CONF->{CAS_URL} || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    return $template->output;
}

sub gen_body {
#    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
#    
#    if ( $USER->is_public ) {
#        $template->param( PAGE_NAME => "$PAGE_TITLE.pl",
#                          LOGIN     => 1 );
#        return $template->output;
#    }

    my $status_code = $FORM->param('status_code');
    
    if ($status_code eq '404') {    
        return '<div style="text-align:center;padding-top:100px"><img src="https://genomevolution.org/coge/picts/gnome.jpg" /><img src="https://genomevolution.org/coge/picts/page_not_found.png" style="vertical-align:top;margin-top:30px;margin-left:-70px;" /></div>';
    }
    else { # 500 and all others
        return '<div style="text-align:center;padding-top:100px"><img src="https://genomevolution.org/coge/picts/gnome.jpg" /><img src="https://genomevolution.org/coge/picts/server_error.png" style="vertical-align:top;margin-top:30px;margin-left:-70px;" /></div>';
    }
}
