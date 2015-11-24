#! /usr/bin/perl -w

use strict;
use CGI;

use CoGeX;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON::XS;

#use DBIx::Class::ResultClass::HashRefInflator;

use vars qw($P $PAGE_TITLE $USER $coge %FUNCTION $FORM $MAX_RESULTS $LINK);

$PAGE_TITLE = 'History';
$MAX_RESULTS = 100;

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

%FUNCTION = (
    get_history_for_user => \&get_history_for_user,
    toggle_star          => \&toggle_star,
    update_comment       => \&update_comment,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    #$template->param( TITLE      => qq{User History} );
    $template->param( PAGE_TITLE => $PAGE_TITLE,
                      TITLE      => 'History',
    				  PAGE_LINK  => $LINK,
    				  HOME       => $P->{SERVER},
                      HELP       => 'History',
                      WIKI_URL   => $P->{WIKI_URL} || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_body {
    if ( $USER->user_name eq 'public' ) {
        my $template =
          HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param( PAGE_NAME => "$PAGE_TITLE" );
        $template->param( LOGIN     => 1 );
        return $template->output;
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( PAGE_NAME  => $PAGE_TITLE . '.pl' );
    $template->param( MAIN       => 1 );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
    $template->param( OPTION_ALL => 1 ) if $USER->is_admin;
    $template->param( USER_NAME  => $USER->name );

   #	$template->param( HISTORY_CONTENTS => get_history_for_user(html_only=>1) );

    return $template->output;
}

sub get_history_for_user {
    my %opts       = @_;
    my $time_range = $opts{time_range};    # in hours
    $time_range = 24 if ( not defined $time_range or $time_range !~ /[-\d]/ );
    my $include_page_accesses = $opts{include_pages_accesses};

    my %users = map { $_->user_id => $_->name } $coge->resultset('User')->all;
    my @entries;
    if ( $USER->is_admin ) {
        if ( $time_range == 0 ) {
            @entries = $coge->resultset('Log')->search(

                #{ description => { 'not like' => 'page access' } },
                { type     => { '!='  => 0 } },
                { order_by => { -desc => 'time' } }
            );
        }
    }
    else {
        if ( $time_range == 0 or $time_range == -3 ) {
            @entries = $coge->resultset('Log')->search(
                {
                    user_id => $USER->id,

                    #{ description => { 'not like' => 'page access' } }
                    type => { '!=' => 0 }
                },
                { order_by => { -desc => 'time' } }
            );
        }
    }

    my @items;
    foreach (@entries) {
        push @items,
          {
            id          => $_->id,
            starred     => ( $_->status != 0 ),
            date_time   => $_->time,
            user        => ( $_->user_id ? $users{ $_->user_id } : 'public' ),
            page        => $_->page,
            description => $_->description,
            link        => ( $_->link ? $_->link : '' ),
            comment     => $_->comment
          };
    }

    # print STDERR "items: " . @items . "\n";
    return encode_json( \@items );
}

sub toggle_star {
    my %opts   = @_;
    my $log_id = $opts{log_id};

    my $entry = $coge->resultset('Log')->find($log_id);
    return '' unless $entry;

    my $status = $entry->status;
    $entry->status( not $status );
    $entry->update();

    return not $status;
}

sub update_comment {
    my %opts    = @_;
    my $log_id  = $opts{log_id};
    my $comment = $opts{comment};

    # print STDERR "udpate_comment: $log_id $comment\n";

    my $entry = $coge->resultset('Log')->find($log_id);
    return unless $entry;

    $entry->comment($comment);
    $entry->update();
}
