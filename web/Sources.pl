#! /usr/bin/perl -w

use strict;
use CGI;
use JSON::XS;
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
  $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE
  $coge $cogeweb %FUNCTION $COOKIE_NAME $FORM $URL
  $COGEDIR $TEMPDIR $TEMPURL);
  
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

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
$TEMPDIR     = $P->{TEMPDIR} . "Sources/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "Sources/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url()) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name => $COOKIE_NAME, coge => $coge) unless $USER;

%FUNCTION = (
	gen_html           => \&gen_html,
	create_source      => \&create_source,
	delete_source      => \&delete_source,
	get_sources        => \&get_sources,
	edit_source_info   => \&edit_source_info,
	update_source_info => \&update_source_info
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

sub create_source {
	my %opts = @_;
	return 0 unless $USER->is_admin;
	return "No specified name!" unless $opts{name};
	
	# Check if one already exists with same name
	my $source = $coge->resultset('DataSource')->find( { name => $opts{name} } );
	return "A data source with the same name already exists" if ($source);
	
	my $link = $opts{link};
	$link =~ s/^\s+//;
	$link = 'http://' . $link if (not $link =~ /^(\w+)\:\/\//);
	
	# Create the new data source
	$coge->resultset('DataSource')->create( 
	  { name => $opts{name},
		description => $opts{desc},
		link => $link
	  } );

	return 1;
}

sub delete_source {
	my %opts = @_;
	return 0 unless $USER->is_admin;
	my $dsid = $opts{dsid};
	return "Must have valid data source id\n" unless ($dsid);
	
	# Delete the data source
	my $ds = $coge->resultset('DataSource')->find($dsid);
	$ds->delete;

	return 1;
}

sub get_sources {
	#my %opts = @_;

	my @sources;
	foreach my $source ( sort {$a->name cmp $b->name } $coge->resultset('DataSource')->all() ) {
		push @sources, 
		  { NAME  => $source->name,
			ID  => $source->id,
			DESC  => $source->description,
			LINK  => $source->link,
			BUTTONS => $USER->is_admin,
			EDIT_BUTTON => "<span class='link ui-icon ui-icon-gear' onclick=\"edit_source_info({dsid: '" . $source->id . "'});\"></span>",
			DELETE_BUTTON => "<span class='link ui-icon ui-icon-trash' onclick=\"delete_source({dsid: '" . $source->id . "'});\"></span>"
		  };
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'Sources.tmpl' );
	$template->param( SOURCE_TABLE => 1 );
	$template->param( SOURCE_LOOP  => \@sources );
	$template->param( BUTTONS      => $USER->is_admin );

	return $template->output;
}

sub edit_source_info {
	my %opts = @_;
	return 0 unless $USER->is_admin;
	my $dsid  = $opts{dsid};
	return 0 unless $dsid;

	my $ds = $coge->resultset('DataSource')->find($dsid);
	my $desc = ( $ds->description ? $ds->description : '' );

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'Sources.tmpl' );
	$template->param( EDIT_SOURCE_INFO => 1 );
	$template->param( DSID             => $dsid );
	$template->param( NAME             => $ds->name );
	$template->param( DESC             => $desc );
	$template->param( LINK             => $ds->link );

	my %data;
	$data{title} = 'Edit Source Info';
	$data{name}   = $ds->name;
	$data{desc}   = $desc;
	$data{link}   = $ds->link;
	$data{output} = $template->output;

	return encode_json( \%data );
}

sub update_source_info {
	my %opts = @_;
	return 0 unless $USER->is_admin;
	my $dsid  = $opts{dsid};
	return 0 unless $dsid;
	my $name = $opts{name};
	return 0 unless $name;
	my $desc = $opts{desc};
	my $link = $opts{link};

	my $ds = $coge->resultset('DataSource')->find($dsid);
	$ds->name($name);
	$ds->description($desc) if $desc;
	$ds->link($link) if ($link);
	$ds->update;

	return 1;
}

sub gen_html {
	my $html;
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP => '/wiki/index.php?title=Lists' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( TITLE      => qq{} );
	$template->param( PAGE_TITLE => qq{Sources} );
	$template->param( LOGO_PNG   => "SourceView-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
	$template->param( BOX_NAME   => " Data Sources:" );
	$template->param( ADJUST_BOX => 1 );
	$html .= $template->output;
}

sub gen_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'Sources.tmpl' );
	$template->param( PAGE_NAME => $FORM->url );
	$template->param( MAIN      => 1 );
	my ($source_info) = get_sources();
	$template->param( SOURCE_INFO => $source_info );
	$template->param( BUTTONS => $USER->is_admin );
	$template->param( ADMIN_AREA => $USER->is_admin );
	return $template->output;
}

