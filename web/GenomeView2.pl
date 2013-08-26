#! /usr/bin/perl -w

use strict;
use CGI;
use HTML::Template;
use CoGeX;
use CoGe::Accessory::Web;
no warnings 'redefine';

use vars
  qw($P $PAGE_TITLE $USER $coge %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS);

$PAGE_TITLE = 'GenomeView2';

$FORM = new CGI;

( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
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

    #	if ($USER->user_name eq 'public') {
    #		$template->param( LOGIN => 1 );
    #		return $template->output;
    #	}

    $template->param( MAIN => 1 );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    #	$template->param( TOC => get_toc(),
    #					  CONTENTS => get_contents(html_only => 1),
    #					  ROLES => get_roles('reader'),
    #					  NOTEBOOK_TYPES => get_notebook_types('mixed') );

    return $template->output;
}

sub get_toc {    # table of contents
    my @rows;
    push @rows,
      {
        TOC_ITEM_ID       => $ITEM_TYPE{mine},
        TOC_ITEM_INFO     => 'My Stuff',
        TOC_ITEM_CHILDREN => 3
      };
    push @rows,
      {
        TOC_ITEM_ID   => $ITEM_TYPE{notebook},
        TOC_ITEM_INFO => 'Notebooks',
        TOC_ITEM_ICON =>
          '<img src="picts/notebook-icon.png" width="15" height="15"/>',
        TOC_ITEM_INDENT => 20
      };
    push @rows,
      {
        TOC_ITEM_ID   => $ITEM_TYPE{genome},
        TOC_ITEM_INFO => 'Genomes',
        TOC_ITEM_ICON =>
          '<img src="picts/dna-icon.png" width="15" height="15"/>',
        TOC_ITEM_INDENT => 20
      };
    push @rows,
      {
        TOC_ITEM_ID   => $ITEM_TYPE{experiment},
        TOC_ITEM_INFO => 'Experiments',
        TOC_ITEM_ICON =>
          '<img src="picts/testtube-icon.png" width="15" height="15"/>',
        TOC_ITEM_INDENT => 20
      };

# push @rows, { TOC_ITEM_ID 	=> $ITEM_TYPE{group},
# 			  TOC_ITEM_INFO 	=> 'Groups',
# 			  TOC_ITEM_ICON 	=> '<img src="picts/group-icon.png" width="15" height="15"/>' };
    push @rows,
      {
        TOC_ITEM_ID   => $ITEM_TYPE{shared},
        TOC_ITEM_INFO => 'Shared with me'
      };
    push @rows,
      {
        TOC_ITEM_ID       => $ITEM_TYPE{activity},
        TOC_ITEM_INFO     => 'Activity',
        TOC_ITEM_CHILDREN => 1
      };
    push @rows,
      {
        TOC_ITEM_ID     => $ITEM_TYPE{activity_viz},
        TOC_ITEM_INFO   => 'Graph',
        TOC_ITEM_INDENT => 20
      };
    push @rows,
      {
        TOC_ITEM_ID   => $ITEM_TYPE{trash},
        TOC_ITEM_INFO => 'Trash'
      };

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        DO_TOC        => 1,
        TOC_ITEM_LOOP => \@rows
    );
    return $template->output;
}

# FIXME these comparison routines are duplicated elsewhere
sub genomecmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->organism->name cmp $b->organism->name
      || versioncmp( $b->version, $a->version )
      || $a->type->id <=> $b->type->id
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

sub experimentcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    versioncmp( $b->version, $a->version )
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

sub listcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    $a->name cmp $b->name;
}

