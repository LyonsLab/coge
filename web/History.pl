#! /usr/bin/perl -w

use strict;
use CGI;

use JSON::XS;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use URI::Escape;
use Data::Dumper;
use File::Path;

no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME 
			$TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION $COOKIE_NAME 
			$FORM $URL $COGEDIR $TEMPDIR $TEMPURL);
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$PAGE_NAME = 'Groups.pl';

$FORM = new CGI;

$DBNAME  = $P->{DBNAME};
$DBHOST  = $P->{DBHOST};
$DBPORT  = $P->{DBPORT};
$DBUSER  = $P->{DBUSER};
$DBPASS  = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL         = $P->{URL};
$COGEDIR     = $P->{COGEDIR};
$TEMPDIR     = $P->{TEMPDIR} . "History/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "History/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

%FUNCTION = (
	gen_html => \&gen_html,
	get_history_for_user => \&get_history_for_user,
);

dispatch();

sub dispatch {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {
		#my %args = $FORM->Vars;
		#print STDERR Dumper \%args;
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else {
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
	else {
		print $FORM->header, gen_html();
	}
}

sub gen_html {
	my $html;
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP => '/wiki/index.php?title=History' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( TITLE      => qq{User History} );
	$template->param( PAGE_TITLE => qq{History} );
	$template->param( LOGO_PNG   => "History-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
#	$name .= $name =~ /s$/ ? "'" : "'s";
#	$template->param( BOX_NAME => $name . " History" );
	$html .= $template->output;
}

sub gen_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'History.tmpl' );
	$template->param( PAGE_NAME => $FORM->url );
	$template->param( MAIN      => 1 );
	$template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
	$template->param( HISTORY_CONTENTS => get_history_for_user() );
	
	return $template->output;
}

sub get_history_for_user
{
	my %opts = @_;
	my $time_range = $opts{time_range};
	$time_range = 1 if (not defined $time_range);
	
	my @entries;
	if ( $USER->is_admin ) {
		if (not $time_range) {
			@entries = $coge->resultset('Log')->all;
		}
		else {
			@entries = $coge->resultset('Log')->search_literal( 'time >= DATE_SUB(NOW(), INTERVAL ? HOUR)', ($time_range) );
		}
	}
	else {
		if (not $time_range) {
			@entries = $coge->resultset('Log')->search( { user_id => $USER->id } );
		}
		else {
			@entries = $coge->resultset('Log')->search_literal( 'user_id = ? AND time >= DATE_SUB(NOW(), INTERVAL ? HOUR)', ($USER->id, $time_range) );
		}
	}

	my @rows;	
	foreach my $entry (reverse @entries) {
		my %row;
		$row{TIME} = $entry->time;
		$row{USER} = ($entry->user ? $entry->user->name : '');
		$row{PAGE} = $entry->page;
		$row{DESCRIPTION} = $entry->description;
		$row{LINK} = '<a href="' . $entry->link . '" target="_blank">' . $entry->link . '</a>' if ($entry->link);
		push @rows, \%row;
	}
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'History.tmpl' );
	$template->param( HISTORY => 1 );
	$template->param( HISTORY_LOOP => \@rows );

	return $template->output;
}
