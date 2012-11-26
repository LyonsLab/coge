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
	gen_html			=> \&gen_html,
	get_logs			=> \&get_logs,
	upload_image_file	=> \&upload_image_file,
);

# debug for fileupload:
# print STDERR $ENV{'REQUEST_METHOD'} . "\n" . $FORM->url . "\n" . Dumper($FORM->Vars) . "\n";	# debug
# print "data begin\n" . $FORM->param('POSTDATA') . "\ndata end\n" if ($FORM->param('POSTDATA'));

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
	return 'Access Denied' if ($USER->user_name eq 'public');

	# Other user specified as param, only allow access if collaborator
	my $uid = $FORM->param('uid');
	my $user = $USER;
	if ($uid) {
		return '' if (not $USER->has_collaborator($uid));
		$user = $coge->resultset('User')->find($uid);
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( PAGE_NAME        => "$PAGE_TITLE.pl" );
	$template->param( MAIN             => 1 );
	$template->param( ADMIN_AREA       => 1 ) if $USER->is_admin;

	my $is_user = ($user->id == $USER->id);
	$template->param( IS_USER => $is_user );
	$template->param( USER_NAME => $user->user_name );
	$template->param( FULL_NAME => $user->display_name );
	$template->param( DESCRIPTION => $user->description );
	$template->param( EMAIL => $user->email );
	$template->param( USER_IMAGE => ($user->image_id ? 'image.pl?id=' . $user->image_id : ($is_user ? 'picts/smiley_default.png' : 'picts/smiley_default2.png') ) );

	$template->param( LOGS => get_logs() );
	$template->param( GROUPS => get_groups($user) );
	$template->param( COLLABORATORS => get_collaborators($user) );
	$template->param( LISTS => get_lists($user) );
	$template->param( GENOMES => get_genomes($user) );
	$template->param( EXPERIMENTS => get_experiments($user) );

	return $template->output;
}

sub get_groups {
	my $user = shift;
	my $is_user = ($user->id == $USER->id);
	my @groups = $user->groups;
	
	my @rows;
	foreach my $group (sort {$a->name cmp $b->name} @groups) {
		next if (!$is_user && !$group->has_member($USER));

		my $id = $group->id;
		my %row;
		$row{GROUP_NAME} = qq{<span class="link" onclick='window.open("GroupView.pl?ugid=$id")'>} . $group->name . "</span>";
		
		my $role = $group->role->name;
#		$role .= ": ".$group->role->description if $group->role->description;
		$row{GROUP_ROLE} = $role;

		$row{GROUP_DESC} = $group->description if $group->description;

		push @rows, \%row;
	}
	push @rows, {GROUP_NAME => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( GROUP_TABLE => 1 );
	$template->param( GROUP_LOOP => \@rows );
	return $template->output;
}

sub get_collaborators {
	my $user = shift;
	my $is_user = ($user->id == $USER->id);
	my @users = $user->collaborators;
	
	my @rows;
	foreach my $u (sort {$a->display_name cmp $b->display_name} @users) {
		next if (!$is_user && !$u->has_collaborator($USER));

		my $id = $u->id;
		my %row;
		$row{COLLABORATOR_INFO} = qq{<span class="link" onclick='window.open("User.pl?uid=$id")'>} . $u->display_name . "</span>";
		
		push @rows, \%row;
	}
	push @rows, {COLLABORATOR_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( COLLABORATOR_TABLE => 1 );
	$template->param( COLLABORATOR_LOOP => \@rows );
	return $template->output;
}

sub get_lists {
	my $user = shift;
	my $is_user = ($user->id == $USER->id);
	my @lists = $user->lists;

	my @rows;
	foreach my $list (sort listcmp @lists) {
		next if ($list->restricted && !$is_user && !$user->has_access(list => $list));

		my $id = $list->id;
		my %row;
		$row{LIST_INFO} = qq{<span class="link" onclick='window.open("NotebookView.pl?lid=$id")'>} . $list->info . "</span>";
		
		push @rows, \%row;
	}
	push @rows, {LIST_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);


	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( LIST_TABLE => 1 );
	$template->param( LIST_LOOP => \@rows );
	return $template->output;
}

sub get_genomes {
	my $user = shift;
	my $is_user = ($user->id == $USER->id);
	my @genomes = $user->owner_group->owner_list->genomes;

	my @rows;
	foreach my $genome (sort genomecmp @genomes) {
		next if ($genome->restricted && !$is_user && !$user->has_access(dsg => $genome));

		my $id = $genome->id;
		my %row;
		$row{GENOME_INFO} = qq{<span class="link" onclick='window.open("OrganismView.pl?dsgid=$id")'>} . $genome->info . "</span>";
		
		push @rows, \%row;
	}
	push @rows, {GENOME_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( GENOME_TABLE => 1 );
	$template->param( GENOME_LOOP => \@rows );
	return $template->output;
}

sub get_experiments {
	my $user = shift;
	my $is_user = ($user->id == $USER->id);
	my @experiments = $user->owner_group->owner_list->experiments;

	my @rows;
	foreach my $exp (sort experimentcmp @experiments) {
		next if ($exp->restricted && !$is_user && !$user->has_access(experiment => $exp));

		my $id = $exp->id;
		my %row;
		$row{EXPERIMENT_INFO} = qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=$id")'>} . $exp->info . "</span>";
		
		push @rows, \%row;
	}
	push @rows, {EXPERIMENT_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( EXPERIMENT_TABLE => 1 );
	$template->param( EXPERIMENT_LOOP => \@rows );
	return $template->output;
}

sub get_logs {
	return if ($USER->user_name eq "public");

	my %opts = @_;
	my $type = $opts{type};

	my @logs;
	if (!$type or $type eq 'recent') {
		@logs = $coge->resultset('Log')->search( { user_id => $USER->id }, { order_by => { -desc => 'time' } } ); # $user->logs;
		#my @logs = reverse $coge->resultset('Log')->search_literal( 'user_id = ' . $user->id . ' AND time >= DATE_SUB(NOW(), INTERVAL 1 HOUR)' );
	}
	else {
		@logs = $coge->resultset('Log')->search( { user_id => $USER->id, status => 1 }, { order_by => { -desc => 'time' } } );
	}

	my @rows;
	foreach (splice(@logs, 0, 100)) {
		push @rows, { LOG_TIME => $_->time,
					  LOG_PAGE => $_->page,
					  LOG_DESC => $_->description,
					  LOG_LINK => $_->link
					};
	}
	return if (not @rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( LOG_TABLE => 1 );
	$template->param( LOG_LOOP => \@rows );
	return $template->output;
}

sub upload_image_file {
	return if ($USER->user_name eq "public");

	my %opts = @_;
	my $image_filename = '' . $FORM->param('input_upload_file');
	my $fh = $FORM->upload('input_upload_file');
	return if (-s $fh > 2*1024*1024); # limit to 2MB

	#TODO delete old image

	# Create the image
	my $image;
	if ($fh) {
		#print STDERR "$image_filename size=" . (-s $fh) . "\n";
		read($fh, my $contents, -s $fh);
		$image = $coge->resultset('Image')->create(
		  {	filename => $image_filename,
			image => $contents
		  }
		);
		return unless $image;

		# Link to user
		$USER->image_id($image->id);
		$USER->update;
		return encode_json({ link => 'image.pl?id=' . $image->id });
	}
	
	return;
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