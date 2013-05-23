#! /usr/bin/perl -w

use strict;
use CGI;

use JSON::XS;
use HTML::Template;
use Sort::Versions;
use List::Util qw(first);
use DBIxProfiler;
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;
use File::Path;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;
use Benchmark;
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE
  $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION
  $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL %ITEM_TYPE
  $MAX_SEARCH_RESULTS);
print STDERR "$ENV{HOME}\n";
$P = CoGe::Accessory::Web::get_defaults("$ENV{HOME}/coge.conf");

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_TITLE = 'GenomeView2';

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
$TEMPDIR     = $P->{TEMPDIR} . "$PAGE_TITLE/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "PAGE_TITLE/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;
my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
$link = CoGe::Accessory::Web::get_tiny_link( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE", url => $link, disable_logging => 1 );

%FUNCTION = (
);

my $node_types = CoGeX::node_types();

%ITEM_TYPE = ( # content/toc types
	all 			=> 100,
	mine			=> 101,
	shared 			=> 102,
	activity 		=> 103,
	trash 			=> 104,
	activity_viz	=> 105,
	user 			=> $node_types->{user},
	group 			=> $node_types->{group},
	notebook 		=> $node_types->{list},
	genome 			=> $node_types->{genome},
	experiment 		=> $node_types->{experiment}
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

sub gen_html {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP       => "/wiki/index.php?title=$PAGE_TITLE",
					  USER       => ($USER->user_name eq "public" ? '' : $USER->display_name),
					  PAGE_TITLE => 'Genome Viewer',
					  LOGO_PNG   => "$PAGE_TITLE-logo.png",
					  DATE       => $DATE,
					  ADJUST_BOX => 1,
					  BODY       => gen_body() );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";

	return $template->output;
}

sub gen_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
	
#	if ($USER->user_name eq 'public') {
#		$template->param( LOGIN => 1 );
#		return $template->output;
#	}

	$template->param( MAIN       => 1 );
	$template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
#	$template->param( TOC => get_toc(),
#					  CONTENTS => get_contents(html_only => 1),
#					  ROLES => get_roles('reader'),
#					  NOTEBOOK_TYPES => get_notebook_types('mixed') );

	return $template->output;
}

sub get_toc { # table of contents
	my @rows;
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{mine}, 
				  TOC_ITEM_INFO 	=> 'My Stuff',
				  TOC_ITEM_CHILDREN => 3 };
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{notebook}, 
				  TOC_ITEM_INFO		=> 'Notebooks', 
				  TOC_ITEM_ICON 	=> '<img src="picts/notebook-icon.png" width="15" height="15"/>', 
				  TOC_ITEM_INDENT 	=> 20 };
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{genome}, 
				  TOC_ITEM_INFO 	=> 'Genomes', 
				  TOC_ITEM_ICON 	=> '<img src="picts/dna-icon.png" width="15" height="15"/>', 
				  TOC_ITEM_INDENT 	=> 20 };
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{experiment}, 
				  TOC_ITEM_INFO 	=> 'Experiments', 
				  TOC_ITEM_ICON 	=> '<img src="picts/testtube-icon.png" width="15" height="15"/>', 
				  TOC_ITEM_INDENT 	=> 20 };
	# push @rows, { TOC_ITEM_ID 	=> $ITEM_TYPE{group}, 
	# 			  TOC_ITEM_INFO 	=> 'Groups', 
	# 			  TOC_ITEM_ICON 	=> '<img src="picts/group-icon.png" width="15" height="15"/>' };
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{shared}, 
	 			  TOC_ITEM_INFO 	=> 'Shared with me' };
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{activity},
				  TOC_ITEM_INFO 	=> 'Activity',
				  TOC_ITEM_CHILDREN => 1 };
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{activity_viz}, 
				  TOC_ITEM_INFO 	=> 'Graph', 
				  TOC_ITEM_INDENT 	=> 20 };				  
	push @rows, { TOC_ITEM_ID 		=> $ITEM_TYPE{trash},
				  TOC_ITEM_INFO 	=> 'Trash' };
	

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_TOC => 1,
					  TOC_ITEM_LOOP => \@rows );
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

