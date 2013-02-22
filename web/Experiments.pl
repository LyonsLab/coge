#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use DBI;

use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON::XS;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use DBIxProfiler;
use File::Path;
use Sort::Versions;
no warnings 'redefine';

use vars qw(
	$P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE 
	$TEMPDIR $USER $DATE $COGEDIR $coge $FORM $URL $TEMPURL $COOKIE_NAME 
	%FUNCTION
);

$P         = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_TITLE = 'Experiments';

$TEMPDIR = $P->{TEMPDIR} . $PAGE_TITLE . '/';
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . $PAGE_TITLE . '/';

$FORM = new CGI;

$DBNAME      = $P->{DBNAME};
$DBHOST      = $P->{DBHOST};
$DBPORT      = $P->{DBPORT};
$DBUSER      = $P->{DBUSER};
$DBPASS      = $P->{DBPASS};
$connstr     = "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge        = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

#$SIG{'__WARN__'} = sub { };    # silence warnings

%FUNCTION = (
	generate_html				=> \&generate_html,
	delete_experiment			=> \&delete_experiment,
	get_experiments_for_user	=> \&get_experiments_for_user,
);

if ( $FORM->param('jquery_ajax') ) {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else {
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
}
else {
	print $FORM->header, "\n", generate_html();
}

sub delete_experiment {
	my %opts = @_;
	my $eid = $opts{eid};
	return "Must have valid experiment id\n" unless ($eid);
	
	# Check permissions
	return unless ($USER->is_admin or $USER->is_owner(experiment => $eid));

	# Delete the experiment and associated connectors & annotations
	#FIXME add some error checking/logging here
	my $experiment = $coge->resultset('Experiment')->find($eid);
	return 0 unless $experiment;
	$experiment->deleted(1);
	$experiment->update;
	
	CoGe::Accessory::Web::log_history( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE", description => 'delete experiment id' . $experiment->id );

	return 1;
}

sub get_experiments_for_user {
	#my %opts = @_;

	my @experiments;
#	if ( $USER->is_admin ) {
#		@experiments = $coge->resultset('Experiment')->all();
#	}
#	else {
		@experiments = $USER->experiments;
#	}

	my @rows;
	foreach my $e (@experiments) {#( sort experimentcmp @experiments ) {
		push @rows, 
		  { NAME  => qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=} . $e->id . qq{")'>} . $e->info . "</span>",
			VERSION  => $e->version,
			DATE =>  $e->date,
			EDIT_BUTTON => 0 ?
				"<span class='link ui-icon ui-icon-locked' onclick=\"alert('This list is locked and cannot be edited.')\"></span>" :
				"<span class='link ui-icon ui-icon-gear' onclick=\"window.open('ExperimentView.pl?eid=" . $e->id . "')\"></span>",
			DELETE_BUTTON => 0 ?
				"<span class='link ui-icon ui-icon-locked' onclick=\"alert('This list is locked and cannot be deleted.')\"></span>" :
				"<span class='link ui-icon ui-icon-trash' onclick=\"dialog_delete_experiment({eid: '" . $e->id . "'});\"></span>"
		  };
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( DO_EXPERIMENT_TABLE => 1 );
	$template->param( EXPERIMENT_LOOP  => \@rows );

	return $template->output;
}

sub generate_html {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( PAGE_TITLE => $PAGE_TITLE );
	$template->param( HELP       => '/wiki/index.php?title=' . $PAGE_TITLE . '.pl' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= ' ' . $USER->last_name if ( $USER->first_name && $USER->last_name );
	$template->param( USER     => $name );
	$template->param( LOGO_PNG => $PAGE_TITLE . "-logo.png" );
	$template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE     => $DATE );
	my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
	$link = CoGe::Accessory::Web::get_tiny_link( url => $link );

	$template->param( BODY => generate_body() );
	$template->param( ADJUST_BOX => 1 );

	return $template->output;
}

sub generate_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( MAIN => 1 );
	$template->param( PAGE_NAME => $PAGE_TITLE );
	$template->param( EXPERIMENT_TABLE  => get_experiments_for_user() );
	
	return $template->output;
}

# FIXME this comparison routine is duplicated elsewhere
sub experimentcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->name cmp $b->name || $b->id cmp $a->id
}
