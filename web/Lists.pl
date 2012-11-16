#! /usr/bin/perl -w

use strict;
use CGI;
#use JSON::XS;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
#use URI::Escape;
use Data::Dumper;
use File::Path;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;

no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS
  $connstr $PAGE_TITLE $TEMPDIR $USER $DATE $BASEFILE
  $coge $cogeweb %FUNCTION $COOKIE_NAME $FORM $URL
  $COGEDIR $TEMPDIR $TEMPURL);
  
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$PAGE_TITLE = 'Lists';

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
($USER) = CoGe::Accessory::Web->login_cas(cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url()) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name => $COOKIE_NAME, coge => $coge) unless $USER;

my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
$link = CoGe::Accessory::Web::get_tiny_link( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE.pl", url => $link );


%FUNCTION = (
	gen_html           => \&gen_html,
	create_list        => \&create_list,
	delete_list        => \&delete_list,
	get_lists_for_user => \&get_lists_for_user,
);

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

sub create_list {
	my %opts = @_;
	return "You need to be a registered user to create a new list!" unless $USER->id;
	return "No specified name!" unless $opts{name};
	return "No specified type!" unless $opts{typeid};
	
	# Get owner user group for the new list
	my $ug = $USER->owner_group;
	return "Internal error!  Please email $P->{SUPPORT_EMAIL} with a problem description.\n" unless $ug;
	
	# Create the new list
	my $list = $coge->resultset('List')->create( 
	  { name => $opts{name},
		description => $opts{desc},
		list_type_id => $opts{typeid},
		user_group_id => $ug->id,
		restricted => 1
	  } );
	
	CoGe::Accessory::Web::log_history( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE.pl", description => 'create list id' . $list->id );

	return 1;
}

sub delete_list {
	my %opts = @_;
	my $lid = $opts{lid};
	return "Must have valid list id\n" unless ($lid);
	
	# Get owner user group for the new list
	my $ug = $USER->owner_group; 
	
	# Delete the list and associated connectors
	#FIXME add some error checking/logging here
	my $list = $coge->resultset('List')->find($lid);
	return 0 if ($list->locked);
	$list->delete;
	
	CoGe::Accessory::Web::log_history( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE.pl", description => 'delete list id' . $list->id );

	return 1;
}

sub get_lists_for_user {
	#my %opts = @_;
	
	my (@lists, @admin_lists);
	if ( $USER->is_admin ) {
	  @admin_lists = $coge->resultset('List')->all();
	}
	@lists = $USER->lists;
	my %user_list_ids = map{$_->id, 1} @lists; #for admins -- marks lists that they don't own
	my %seen_list_ids; #for admins -- their lists are listed first, then all lists
	my @list_info;
	@lists = sort sort {$a->name cmp $b->name } @lists;
	push @lists, sort sort {$a->name cmp $b->name } @admin_lists;
	foreach my $list (@lists ) {
	  next if $seen_list_ids{$list->id};
	  $seen_list_ids{$list->id}=1;
	  my $name = qq{<span class=link onclick='window.open("ListView.pl?lid=} . $list->id . qq{")'>} . $list->name . "</span>";
	  $name = '&reg; '.$name if $list->restricted;
	  $name = '&alpha; '.$name if !$user_list_ids{$list->id};

	  push @list_info, 
	    { 
	     NAME  => $name,
	     DESC  => $list->description,
	     TYPE  => ( $list->type ? $list->type->name : '' ),
	     ANNO  => join( "<br>", map { $_->type->name . ": " . $_->annotation . ($_->image ? ' (image)' : '') } sort sortAnno $list->annotations ),
	     DATA  => $list->data_summary(),
	     GROUP => $list->group->info_html,
	     RESTRICTED => ($list->restricted ? 'yes' : 'no'),
	     EDIT_BUTTON => $list->locked ?
	     "<span class='link ui-icon ui-icon-locked' onclick=\"alert('This list is locked and cannot be edited.')\"></span>" :
	     "<span class='link ui-icon ui-icon-gear' onclick=\"window.open('ListView.pl?lid=" . $list->id . "')\"></span>",
	     DELETE_BUTTON => $list->locked ?
	     "<span class='link ui-icon ui-icon-locked' onclick=\"alert('This list is locked and cannot be deleted.')\"></span>" :
	     "<span class='link ui-icon ui-icon-trash' onclick=\"dialog_delete_list({lid: '" . $list->id . "'});\"></span>"
	    };
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( LIST_STUFF => 1 );
	$template->param( LIST_LOOP  => \@list_info );
	$template->param( TYPE_LOOP => get_list_types() );

	return $template->output;
}

sub sortAnno {
	$a->type->name cmp $b->type->name || $a->annotation cmp $b->annotation
}

sub get_list_types {
	my @types;
	foreach my $type ( $coge->resultset('ListType')->all() ) {
		next if ($type->name =~ /owner/i); # reserve this type for system-created lists
		my $name = $type->name . ($type->description ? ": " . $type->description : '');
		push @types, { TID => $type->id, NAME => $name };
	}
	return \@types;
}

sub gen_html {
	my $html;
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( TITLE      => qq{} );
	$template->param( PAGE_TITLE => $PAGE_TITLE );
	$template->param( LOGO_PNG   => "$PAGE_TITLE-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
#	$name .= $name =~ /s$/ ? "'" : "'s";
#	$template->param( BOX_NAME   => $name . " Data Lists:" );
	$template->param( ADJUST_BOX => 1 );
	$html .= $template->output;
}

sub gen_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( PAGE_NAME  => "$PAGE_TITLE.pl" );
	$template->param( MAIN       => 1 );
	$template->param( LIST_INFO  => get_lists_for_user() );
	$template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
	return $template->output;
}

