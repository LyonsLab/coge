#! /usr/bin/perl -w

use strict;
use CGI;

use JSON::XS;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use URI::Escape;
use Data::Dumper;
use File::Path;
use DBIx::Class::ResultClass::HashRefInflator;

no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE 
			$TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION $COOKIE_NAME 
			$FORM $URL $COGEDIR $TEMPDIR $TEMPURL $MAX_RESULTS);
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$PAGE_TITLE = 'History';

$MAX_RESULTS = 100;

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
$TEMPDIR     = $P->{TEMPDIR} . "$PAGE_TITLE/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "$PAGE_TITLE/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

%FUNCTION = (
	gen_html 				=> \&gen_html,
	get_history_for_user	=> \&get_history_for_user,
	toggle_star				=> \&toggle_star,
	update_comment			=> \&update_comment,
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
	$template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( TITLE      => qq{User History} );
	$template->param( PAGE_TITLE => $PAGE_TITLE );
	$template->param( LOGO_PNG   => "$PAGE_TITLE-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
#	$name .= $name =~ /s$/ ? "'" : "'s";
#	$template->param( BOX_NAME => $name . " History" );
	$html .= $template->output;
}

sub gen_body {
	if ($USER->user_name eq 'public') {
		my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
		$template->param( PAGE_NAME => "$PAGE_TITLE" );
		$template->param( LOGIN     => 1 );
		return $template->output;
	}
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( PAGE_NAME  => $PAGE_TITLE . '.pl' );
	$template->param( MAIN       => 1 );
	$template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
	$template->param( USER_NAME  => $USER->name );
#	$template->param( HISTORY_CONTENTS => get_history_for_user(html_only=>1) );
	
	return $template->output;
}

sub get_history_for_user
{
	my %opts = @_;
	my $time_range = $opts{time_range}; # in hours
	$time_range = 24 if (not defined $time_range or $time_range !~ /[-\d]/);
	my $last_id = $opts{last_id};
#	print STDERR "get_history_for_user: $time_range $last_id\n";
	
	my %users = map { $_->user_id => $_->name } $coge->resultset('User')->all;
	my @entries;
	if ( $USER->is_admin ) {
		if ($time_range == 0) {
			#$rs = $coge->resultset('Log')->search( undef, { result_class => 'DBIx::Class::ResultClass::HashRefInflator', order_by => { -desc => 'time' } });
			@entries = $coge->resultset('Log')->search( undef, { order_by => { -desc => 'time' } });
		}
	}
	else {
		if ($time_range == 0 or $time_range == -3) {
			@entries = $coge->resultset('Log')->search( { user_id => $USER->id }, { order_by => { -desc => 'time' } } );
		}
	}

	my @items;
	foreach (@entries) {
		push @items, 
			{ id => $_->id,
			  starred => ($_->status != 0),
			  date_time => $_->time, 
			  user => ($_->user_id ? $users{$_->user_id} : 'public'),
			  page => $_->page,
			  description => $_->description,
			  link => ($_->link ? $_->link : ''),
			  comment => $_->comment
			};
	}

	# print STDERR "items: " . @items . "\n";
	return encode_json(\@items);
}

sub toggle_star 
{
	my %opts = @_;
	my $log_id = $opts{log_id};
	
	my $entry = $coge->resultset('Log')->find($log_id);
	return '' unless $entry;
	
	my $status = $entry->status;
	$entry->status(not $status);
	$entry->update();
	
	return not $status;
}

sub update_comment
{
	my %opts = @_;
	my $log_id = $opts{log_id};
	my $comment = $opts{comment};
	# print STDERR "udpate_comment: $log_id $comment\n";

	my $entry = $coge->resultset('Log')->find($log_id);
	return unless $entry;
	
	$entry->comment($comment);
	$entry->update();
}
