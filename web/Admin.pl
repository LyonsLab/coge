#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON qw(encode_json);
no warnings 'redefine';

use vars qw($P $PAGE_NAME $USER $BASEFILE $coge $cogeweb %FUNCTION $FORM $MAX_SEARCH_RESULTS);

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init( cgi => $FORM );

$MAX_SEARCH_RESULTS = 400;

%FUNCTION = (
    search_organisms     => \&search_organisms,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => '/wiki/index.php?title=ADMIN' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER       => $name );
    $template->param( PAGE_TITLE => qq{Admin} );
    $template->param( LOGO_PNG   => "Blank-logo.png" );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    $template->param( ADJUST_BOX => 1 );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'Admin.tmpl' );
    $template->param( MAIN => 1);

    #print STDERR "Can you see this?";
    #print STDERR $template->output;
    return $template->output;
}

sub search_organisms {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @organisms = $coge->resultset("Organism")->search(
        \[
            'name LIKE ? OR description LIKE ?'
            ,    #FIXME security hole: need to check 'restricted'
            [ 'name', $search_term ], [ 'description', $search_term ]
        ]
    );

    # Limit number of results displayed
    if ( @organisms > $MAX_SEARCH_RESULTS ) {
        return encode_json( { timestamp => $timestamp, items => undef } );
    }

    my @results;
    foreach ( sort { $a->name cmp $b->name } @organisms ) {
        push @results, { 'label' => $_->name, 'value' => $_->id };
    }

    return encode_json( { timestamp => $timestamp, items => \@results } );
}

sub search_users {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @users = $coge->resultset("User")->search(
        \[
            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
            [ 'user_name',  $search_term ],
            [ 'first_name', $search_term ],
            [ 'last_name',  $search_term ]
        ]
    );

    # Limit number of results displayed
    # if (@users > $MAX_SEARCH_RESULTS) {
    # 	return encode_json({timestamp => $timestamp, items => undef});
    # }

    return encode_json(
        {
            timestamp => $timestamp,
            items     => [ sort map { $_->user_name } @users ]
        }
    );
}
