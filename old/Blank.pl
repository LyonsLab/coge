#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use HTML::Template;
no warnings 'redefine';

use vars qw($P $PAGE_NAME $USER $BASEFILE $coge $cogeweb %FUNCTION $FORM);

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init( cgi => $FORM );

%FUNCTION = (
    gen_html  => \&gen_html,
    test_ajax => \&test_ajax,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => '/wiki/index.php?title=BLANK' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER       => $name );
    $template->param( TITLE      => qq{Blank Page} );
    $template->param( PAGE_TITLE => qq{Blank Page} );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'Blank.tmpl' );
    $template->param( PAGE_NAME => $FORM->url );

    return $template->output;
}

sub test_ajax {
    my %opts   = @_;
    my $output = "Test worked!\n";
    $output .= "<pre>" . "Args: " . Dumper( \%opts ) . "</pre>";
}
