#! /usr/bin/perl -w

use strict;
use CGI;

#use CGI::Ajax;
use JSON::XS;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use Sort::Versions;
use List::Util qw(first);
use DBIxProfiler;
#use URI::Escape;
use Data::Dumper;
use File::Path;
use File::stat;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;
use Benchmark;
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE
  $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION
  $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL );
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_TITLE = 'User';

$FORM = new CGI;

$DBNAME  = $P->{DBNAME};
$DBHOST  = $P->{DBHOST};
$DBPORT  = $P->{DBPORT};
$DBUSER  = $P->{DBUSER};
$DBPASS  = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL         = $P->{URL};
$COGEDIR     = $P->{COGEDIR};
$TEMPDIR     = $P->{TEMPDIR} . "PAGE_TITLE/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "PAGE_TITLE/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;
my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
$link = CoGe::Accessory::Web::get_tiny_link( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE.pl", url => $link );


%FUNCTION = (
	gen_html	=> \&gen_html,
);

# debug for fileupload:
#print STDERR $ENV{'REQUEST_METHOD'} . "\n" . $FORM->url . "\n" . Dumper($FORM->Vars) . "\n";	# debug
#print "data begin\n" . $FORM->param('POSTDATA') . "\ndata end\n" if ($FORM->param('POSTDATA'));

dispatch();

sub dispatch {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {
		die if not defined $FUNCTION{$fname};
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
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP       => "/wiki/index.php?title=$PAGE_TITLE" );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( TITLE      => 'User Profile' );
	$template->param( PAGE_TITLE => 'User Profile' );
	$template->param( LOGO_PNG   => "$PAGE_TITLE-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
	$template->param( ADJUST_BOX => 1 );

	return $template->output;
}

sub gen_body {
	return "Access Denied" if ($USER->user_name eq "public");

	my $uid = $FORM->param('uid');
	my $user = $USER;
	if ($uid) {
		$user = $coge->resultset('User')->find($uid);
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( PAGE_NAME        => "$PAGE_TITLE.pl" );
	$template->param( MAIN             => 1 );
	$template->param( ADMIN_AREA       => 1 ) if $USER->is_admin;

	$template->param( CAN_EDIT => ($user == $USER) );
	$template->param( USER_NAME => $user->user_name );
	$template->param( FULL_NAME => $user->display_name );
	$template->param( DESCRIPTION => $user->description );
	$template->param( EMAIL => $user->email );
	$template->param( HISTORY_COUNT => $user->history(count=>1) );

	$template->param( GROUPS => get_groups($user) );
	$template->param( COLLABORATORS => get_collaborators($user) );
	$template->param( LISTS => get_lists($user) );
	$template->param( GENOMES => get_genomes($user) );
	$template->param( EXPERIMENTS => get_experiments($user) );

	return $template->output;
}

sub get_groups {
	my $user = shift;
	my @groups = $user->groups;
	
	my @rows;
	foreach my $group (sort {$a->name cmp $b->name} @groups) {
		my $id = $group->id;
				
		my %row;
		$row{GROUP_NAME} = qq{<span class="link" onclick='window.open("GroupView.pl?ugid=$id")'>} . $group->name . "</span>";
		
		my $role = $group->role->name;
#		$role .= ": ".$group->role->description if $group->role->description;
		$row{GROUP_ROLE} = $role;

		$row{GROUP_DESC} = $group->description if $group->description;

		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( GROUP_TABLE => 1 );
	$template->param( GROUP_LOOP => \@rows );
	return $template->output;
}

sub get_collaborators {
	my $user = shift;
	my @users = $user->collaborators;
	
	my @rows;
	foreach my $u (sort {$a->display_name cmp $b->display_name} @users) {
		my $id = $u->id;
				
		my %row;
		$row{COLLABORATOR_INFO} = qq{<span class="link" onclick='window.open("User.pl?uid=$id")'>} . $u->display_name . "</span>";
		
		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( COLLABORATOR_TABLE => 1 );
	$template->param( COLLABORATOR_LOOP => \@rows );
	return $template->output;
}

sub get_lists {
	my $user = shift;
	my @lists = $user->lists;
	
	my @rows;
	foreach my $list (sort listcmp @lists) {
		my $id = $list->id;
				
		my %row;
		$row{LIST_INFO} = qq{<span class="link" onclick='window.open("ListView.pl?lid=$id")'>} . $list->info . "</span>";
		
		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( LIST_TABLE => 1 );
	$template->param( LIST_LOOP => \@rows );
	return $template->output;
}

sub get_genomes {
	my $user = shift;
	my @genomes = $user->owner_group->owner_list->genomes;
	
	my @rows;
	foreach my $genome (sort genomecmp @genomes) {
		my $id = $genome->id;
				
		my %row;
		$row{GENOME_INFO} = qq{<span class="link" onclick='window.open("OrganismView.pl?dsgid=$id")'>} . $genome->info . "</span>";
		
		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( GENOME_TABLE => 1 );
	$template->param( GENOME_LOOP => \@rows );
	return $template->output;
}

sub get_experiments {
	my $user = shift;
	my @experiments = $user->owner_group->owner_list->experiments;
	
	my @rows;
	foreach my $exp (sort experimentcmp @experiments) {
		my $id = $exp->id;
				
		my %row;
		$row{EXPERIMENT_INFO} = qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=$id")'>} . $exp->info . "</span>";
		
		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( EXPERIMENT_TABLE => 1 );
	$template->param( EXPERIMENT_LOOP => \@rows );
	return $template->output;
}

# FIXME these comparison routines are duplicated elsewhere
sub genomecmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	$a->organism->name cmp $b->organism->name || versioncmp($b->version, $a->version) || $a->type->id <=> $b->type->id || $a->name cmp $b->name || $b->id cmp $a->id
}

sub experimentcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->name cmp $b->name || $b->id cmp $a->id
}

sub listcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	$a->name cmp $b->name
}