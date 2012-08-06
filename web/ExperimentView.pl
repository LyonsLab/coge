#! /usr/bin/perl -w

use strict;
use CGI;
use lib '/home/mbomhoff/CoGeX/lib';    # FIXME mdb remove this before commit
use CoGeX;
use DBI;

use Data::Dumper;
use lib '/home/mbomhoff/CoGe/Accessory/lib';    # FIXME mdb remove this before commit
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use DBIxProfiler;
use File::Path;
no warnings 'redefine';

use vars qw( $P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME
  $TEMPDIR $USER $DATE $COGEDIR $coge $FORM $URL
  $TEMPURL $COOKIE_NAME $EXPID );

$P         = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE      = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$PAGE_NAME = "ExperimentView.pl";

$TEMPDIR = $P->{TEMPDIR} . "ExperimentView/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "ExperimentView/";

$FORM = new CGI;

$EXPID = $FORM->param('expid');

$DBNAME  = $P->{DBNAME};
$DBHOST  = $P->{DBHOST};
$DBPORT  = $P->{DBPORT};
$DBUSER  = $P->{DBUSER};
$DBPASS  = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge    = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
	ticket   => $cas_ticket,
	coge     => $coge,
	this_url => $FORM->url()
  )
  if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
	cookie_name => $COOKIE_NAME,
	coge        => $coge
  )
  unless $USER;

$SIG{'__WARN__'} = sub { };    # silence warnings

my %FUNCTION = ( generate_html => \&generate_html,
				 remove_group => \&remove_group,
				 get_groups => \&get_groups );

if ( $FORM->param('jquery_ajax') ) {
	dispatch();
}
else {
	print $FORM->header, "\n", generate_html();
}

sub dispatch {
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

sub remove_group
{
	my %opts  = @_;
	my $ugid  = $opts{ugid};
	my $expid = $opts{expid};
	my $ugec = $coge->resultset('UserGroupExperimentConnector')->find({ user_group_id => $ugid, experiment_id => $expid });
	$ugec->delete();
}

sub generate_html {
	my $html;
	my ($body) = generate_body();
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( PAGE_TITLE => 'ExperimentView' );
	$template->param( HELP       => '/wiki/index.php?title=ExperimentView' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= ' ' . $USER->last_name
	  if ( $USER->first_name && $USER->last_name );
	$template->param( USER     => $name );
	$template->param( LOGO_PNG => "ExperimentView-logo.png" );
	$template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE     => $DATE );
	my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
	$link = CoGe::Accessory::Web::get_tiny_link( url => $link );
	my $box_name = "Experiment List: ";
	my $list_name = $FORM->param('list_name') || $FORM->param('ln');
	$box_name .= " $list_name" if $list_name;
	$box_name .= "<span class=link onclick=window.open('$link');>$link</span>";
	$template->param( BOX_NAME   => $box_name );
	$template->param( BODY       => $body );
	$template->param( ADJUST_BOX => 1 );
	$html .= $template->output;
}

sub get_groups {
	my %opts  = @_;
	my $expid = $opts{expid};
	my $groups;
	
	my $exp = $coge->resultset('Experiment')->find($expid);
	if ($exp) {	
		foreach my $ug ( $exp->user_groups ) {
			$groups .= "<tr>";
			
			$groups .= "<td>" . $ug->name . ( $ug->description ? ': ' . $ug->description : '' ) . "</td>";
			$groups .= "<td>" . $ug->role->name . "</td>";
			$groups .= "<td>" . join( ", ", map { $_->name } $ug->role->permissions ) . "</td>";
	
			my @users;
			foreach my $user ( $ug->users ) {
				push @users, $user->name;
			}
			$groups .= "<td>" . join( ",<br>", @users ) . "</td>";
			my $ugid = $ug->id;
			$groups .= qq{<td><span class='ui-button ui-corner-all ui-button-icon-left'><span class="ui-icon ui-icon-trash" onclick="} . "remove_group({ugid: $ugid, expid: $expid});" . qq{"</span></span></td>};
			
			$groups .= "</tr>";
		}
	}
    return $groups;	
}

sub generate_body {
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentView.tmpl' );
	my $form = $FORM;

	my $exp = $coge->resultset('Experiment')->find($EXPID);
	if ($exp) {
		# Create make private/public button
		my $access_button;
		my $edit_button;
		if (1) {      #$USER->is_owner(dsg=>$dsgid) || $USER->is_admin) {
			if ( $exp->restricted ) {
				$access_button = qq{<span class="ui-button ui-corner-all ui-button-go" onClick="make_experiment_public('$EXPID')">Make Experiment Public</span>};
			}
			else {
				$access_button = qq{<span class="ui-button ui-corner-all ui-button-go" onClick="make_experiment_private('$EXPID')">Make Experiment Private</span>};
			}
			$edit_button = qq{<span class="ui-button ui-corner-all ui-button-go" onClick="edit_experiment_info('$EXPID')">Edit Experiment Info</span>}; 
		}

		# Build sub-table of types
		my %types;
		foreach ( $exp->types ) {
			my $key = $_->name . ( $_->desc ? ': ' . $_->desc : '' );
			$types{$key}++;
		}
		my $type_tbl = '<table class="small"><tbody>';
		foreach ( sort keys %types ) {
			$type_tbl .= "<tr><td>$_</td></tr>";
		}
		$type_tbl .= '</tbody></table>';

		# Build source link
		my $src  = $exp->source;
		my $link = $src->link;
		$link = 'http://' . $link if ( not $link =~ /^http:\/\// );
		my $source = "<span class=link onclick=window.open('$link');>" . $src->name . ': ' . $src->desc . "</span>";

		# Build sub-table of annotations
		%types = ();
		foreach ( $exp->annotations ) {
			my $type      = $_->annotation_type;
			my $group     = $type->group;
			my $groupname = ( $group ? $group->name : '' );
			my $typename  = $type->name;
			push @{ $types{$groupname}{$typename} }, $_->annotation;
		}
		my $annot_tbl;
		foreach my $groupname ( sort keys %types ) {
			my $first1 = 1;
			foreach my $typename ( sort keys %{ $types{$groupname} } ) {
				my $first2 = 1;
				foreach my $annot ( sort @{ $types{$groupname}{$typename} } ) {
					$annot_tbl .= '<tr>';
					$annot_tbl .=
					  '<td>' . ( $first1 ? $groupname : '' ) . '</td>'
					  if ($groupname);
					$annot_tbl .=
					  '<td>' . ( $first2 ? $typename : '' ) . '</td>';
					$annot_tbl .= '<td></td>' if ( not $groupname );
					$annot_tbl .= '<td>' . $annot . '</td></tr>';
					$first1 = $first2 = 0;
				}
			}
		}

		# Build sub-table of user groups
		my $groups = get_groups(expid => $EXPID);

		# Populate the template
		$template->param( EXP_NAME      => $exp->name );
		$template->param( EXP_DESC      => $exp->desc );
		$template->param( EXP_TYPE      => $type_tbl );
		$template->param( EXP_SOURCE    => $source );
		$template->param( EXP_VER       => $exp->version );
		$template->param( ACCESS_BUTTON => $access_button );
		$template->param( EDIT_BUTTON 	=> $edit_button );
		$template->param( EXP_ANNOT     => $annot_tbl );
		$template->param( GROUPS		=> $groups);
	}

	#	if ($table) {
	return $template->output;

	#	}
	#	else {
	#		return "No dataset_group (genome) or experiment ids were specified.";
	#	}
}
