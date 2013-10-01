#! /usr/bin/perl -w

use strict;
use CGI;
use HTML::Template;
use CoGeX;
use CoGe::Accessory::Web;
no warnings 'redefine';

use vars
  qw($P $PAGE_TITLE $USER $coge %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS $LINK);

$PAGE_TITLE = 'GenomeView';

$FORM = new CGI;

( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    ticket     => $FORM->param('ticket') || undef,
    url        => $FORM->url,
    page_title => $PAGE_TITLE
);

%FUNCTION = ();

my $node_types = CoGeX::node_types();

%ITEM_TYPE = (    # content/toc types
    all          => 100,
    mine         => 101,
    shared       => 102,
    activity     => 103,
    trash        => 104,
    activity_viz => 105,
    user         => $node_types->{user},
    group        => $node_types->{group},
    notebook     => $node_types->{list},
    genome       => $node_types->{genome},
    experiment   => $node_types->{experiment}
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param(
        HELP => "/wiki/index.php?title=$PAGE_TITLE",
        USER => ( $USER->user_name eq "public" ? '' : $USER->display_name ),
        PAGE_TITLE => 'Genome Viewer',
        PAGE_LINK  => $LINK,
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

    $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
    $template->param( MAIN => 1 );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}
