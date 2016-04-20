#! /usr/bin/perl -w

use strict;
use CoGe::Accessory::Web;
use CGI;
use HTML::Template;

use vars qw($P $PAGE_TITLE
  $USER $BASEFILE $coge $cogeweb %FUNCTION
  $FORM    );

$PAGE_TITLE = 'SeqType';

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
    $template->param( HELP => '/wiki/index.php?title=SeqType' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(
        USER       => $name,
        PAGE_TITLE => $PAGE_TITLE,
    );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    return $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SeqType.tmpl' );
    $template->param( PAGE_NAME => $FORM->url );
    my @types;
    foreach my $st ( $coge->resultset('GenomicSequenceType')->all() ) {
        my $count = $st->genomes->count();
        push @types,
          {
            ID     => $st->id,
            NAME   => $st->name,
            DESC   => $st->description,
            GCOUNT => $count,
          };
    }
    $template->param( SEQTYPE => \@types );
    return $template->output;
}
