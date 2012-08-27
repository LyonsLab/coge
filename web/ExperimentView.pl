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

use vars qw( $P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME
	$TEMPDIR $USER $DATE $COGEDIR $coge $FORM $URL
	$TEMPURL $COOKIE_NAME %FUNCTION);

$P         = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$PAGE_NAME = "ExperimentView.pl";

$TEMPDIR = $P->{TEMPDIR} . "ExperimentView/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "ExperimentView/";

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

$SIG{'__WARN__'} = sub { };    # silence warnings

%FUNCTION = (
	generate_html           => \&generate_html,
	remove_group            => \&remove_group,
	get_groups              => \&get_groups,
	get_experiment_info     => \&get_experiment_info,
	edit_experiment_info    => \&edit_experiment_info,
	update_experiment_info  => \&update_experiment_info,
	make_experiment_public  => \&make_experiment_public,
	make_experiment_private => \&make_experiment_private,
	search_genomes          => \&search_genomes,
	add_experiment_type     => \&add_experiment_type,
	add_type_to_experiment  => \&add_type_to_experiment,
	get_experiment_types    => \&get_experiment_types,
	get_type_description    => \&get_type_description,
	remove_experiment_type  => \&remove_experiment_type,
	get_experiment_annotations   => \&get_experiment_annotations,
	add_experiment_annotation    => \&add_experiment_annotation,
	add_annotation_to_experiment => \&add_annotation_to_experiment,
	remove_experiment_annotation => \&remove_experiment_annotation,
	search_annotation_types      => \&search_annotation_types,
	get_annotation_type_groups   => \&get_annotation_type_groups,
);

if ( $FORM->param('jquery_ajax') ) {
	dispatch();
}
else {
	print $FORM->header, "\n", generate_html();
}

sub dispatch {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname)
	{
		if ( $args{args} )
		{
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else
		{
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
}

sub remove_group {
	my %opts  = @_;
	my $ugid  = $opts{ugid};
	my $expid = $opts{expid};
	my $ugec  = $coge->resultset('UserGroupExperimentConnector')->find( { user_group_id => $ugid, experiment_id => $expid } );
	$ugec->delete();
}

sub edit_experiment_info {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;

	my $exp = $coge->resultset('Experiment')->find($eid);
	my $desc = ( $exp->description ? $exp->description : '' );

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentView.tmpl' );
	$template->param( EDIT_EXPERIMENT_INFO => 1 );
	$template->param( EID            => $eid );
	$template->param( NAME           => $exp->name );
	$template->param( DESC           => $desc );
	$template->param( GENOME         => $exp->genome->info );
	$template->param( GENOME_ID      => $exp->genome->id );
	$template->param( SOURCE_LOOP    => get_sources($exp->source->id) );
	$template->param( VERSION        => $exp->version );

	my %data;
	$data{title} = 'Edit Experiment Info';
	$data{name}   = $exp->name;
	$data{desc}   = $desc;
	$data{output} = $template->output;

	return encode_json( \%data );
}

sub update_experiment_info {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;
	my $name = $opts{name};
	return 0 unless $name;
	my $desc = $opts{desc};
	my $genome_id = $opts{genome_id};
	my $source_id = $opts{source_id};
	my $version = $opts{version};

	my $exp = $coge->resultset('Experiment')->find($eid);
	$exp->name($name);
	$exp->description($desc) if $desc;
	$exp->version($version);
	$exp->genome_id($genome_id);
	$exp->data_source_id($source_id);
	$exp->update;

	return 1;
}

sub make_experiment_public {
	my %opts  = @_;
	my $eid = $opts{eid};
	return "No EID specified" unless $eid;
	#return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $exp = $coge->resultset('Experiment')->find($eid);
	$exp->restricted(0);
	$exp->update;

	return 1;
}

sub make_experiment_private {
	my %opts  = @_;
	my $eid = $opts{eid};
	return "No EID specified" unless $eid;
	#return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	my $exp = $coge->resultset('Experiment')->find($eid);
	$exp->restricted(1);
	$exp->update;

	return 1;
}

sub genomecmp { # FIXME mdb 8/24/12 - redundant declaration in ListView.pl
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->type->id <=> $b->type->id || $a->name cmp $b->name || $b->id cmp $a->id
}

sub search_genomes {
	my %opts = @_;
	my $search_term = $opts{search_term};
#	print STDERR "$search_term\n";
	return '' unless $search_term and length $search_term > 0;
	
	$search_term = '%'.$search_term.'%';
	
	my @organisms = $coge->resultset("Organism")->search(
		\[ 'name LIKE ? OR description LIKE ?', 
		['name', $search_term ], ['description', $search_term] ]);

	my @genomes = $coge->resultset("Genome")->search(
		\[ 'restricted=? AND (name LIKE ? OR description LIKE ?)', 
		['restricted', 0], ['name', $search_term ], ['description', $search_term] ]);

	# Combine genomes with organism search but prevent duplicates
	# FIXME mdb 8/21/12 - could be done more elegantly
	my %gid;
	map { $gid{$_->id}++ } @genomes;
	foreach (@organisms) {
		foreach ($_->genomes) {
			push @genomes, $_ if (not defined $gid{$_->id});
		}
	}
	
	my @results;
	foreach (sort genomecmp @genomes) {
		push @results, { 'label' => $_->info, 'value' => $_->id }; 
	}
	
	return encode_json( \@results );
}

sub get_sources {
	my $current_source_id = shift;
	
	my @sources;
	foreach my $source ( $coge->resultset('DataSource')->all() ) {
		my $name = $source->name . ($source->description ? ": " . $source->description : '');
		my $selected = '';
		$selected = 'selected="selected"' if ($source->id == $current_source_id);
		push @sources, { SID => $source->id, NAME => $name, SOURCE_SELECTED => $selected };
	}
	return \@sources;
}

sub add_experiment_type {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;

	my $exp = $coge->resultset('Experiment')->find($eid);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentView.tmpl' );
	$template->param( ADD_EXPERIMENT_TYPE => 1 );
	$template->param( EID => $eid );

	my %data;
	$data{title} = 'Add Type';
	$data{output} = $template->output;

	return encode_json( \%data );
}

sub add_type_to_experiment {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;
	my $name = $opts{name};
	return 0 unless $name;
	my $description = $opts{description};
	
	my $type = $coge->resultset('ExperimentType')->find( { name => $name } );
	
	if ($type) {
		# If type exists, check if already assigned to this experiment
		foreach ($type->experiment_type_connectors) {
			if ($_->experiment_id == $eid) {
				if ($type->description ne $description) {
					$type->description($description);
					$type->update;
				}
				return 1;
			}	
		}
	}
	else {
		# Create type if it doesn't already exist
		$type = $coge->resultset('ExperimentType')->create( { name => $name, description => $description } );
	}
	
	# Create connection
	$coge->resultset('ExperimentTypeConnector')->create( 
		{ experiment_id => $eid, 
		  experiment_type_id => $type->id
		} );
	
	return 1;
}

sub get_experiment_types {
	#my %opts = @_;
	
	my %unique;
	
	my $rs = $coge->resultset('ExperimentType');
	while (my $et = $rs->next) {
		$unique{$et->name}++;
	}
	
	return encode_json( [ sort keys %unique ] );
}

sub get_type_description {
	my %opts = @_;
	my $name = $opts{name};
	return unless ($name);
	
	my $type = $coge->resultset('ExperimentType')->find( { name => $name } );
	return $type->description if ($type);
}

sub linkify {
	my ($link, $desc) = @_;
	return "<span class='link' onclick=\"window.open('$link')\">" . $desc . "</span>";
}

sub remove_experiment_type {
	my %opts  = @_;
	my $eid = $opts{eid};
	return "No experiment ID specified" unless $eid;
	my $etid = $opts{etid};
	return "No experiment type ID specified" unless $etid;	
	#return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $etc = $coge->resultset('ExperimentTypeConnector')->find( { experiment_id => $eid, experiment_type_id => $etid } );
	$etc->delete();
	
	return 1;
}

sub get_experiment_annotations {
	my %opts = @_;
	my $eid  = $opts{eid};
	return "Must have valid experiment id\n" unless ($eid);
	return "Access denied\n" unless ($USER->has_access(experiment=>$eid));
	
	my ($exp) = $coge->resultset('Experiment')->find($eid);
	return unless $exp;
	
	my $user_can_edit = ($USER->is_admin || $USER->is_owner_editor(experiment => $eid));
	
	my $anno_obj = new CoGe::Accessory::Annotation( Type => 'anno' );
	$anno_obj->Val_delimit("\n");
	$anno_obj->Add_type(0);
	$anno_obj->String_end("\n");
	
	my %groups;
	foreach my $a ( $exp->annotations ) {
		my $group = (defined $a->type->group ? $a->type->group->name . ':' . $a->type->name : $a->type->name);
		push @{$groups{$group}}, $a;
	}
	
	foreach my $group (sort keys %groups) {
		foreach my $a ( @{$groups{$group}} ) { #( $list->annotations ) {
			my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . $group . "</span>" );
			$anno_type->Type_delimit(": <td class=\"data5\">");
			$anno_type->Val_delimit("<br>");
			my $a_info = ($a->link ? linkify($a->link, $a->info) : $a->info);
			$a_info .= ($user_can_edit ? "<span onClick=\"remove_experiment_annotation({eid: '$eid', eaid: '" . $a->id . "'});\" class=\"link ui-icon ui-icon-trash\"></span>" : '');
			$anno_type->add_Annot($a_info);
			$anno_obj->add_Annot($anno_type);
		}
	}
	my $html = "<table cellpadding=0 class='ui-widget-content ui-corner-all small'>" . $anno_obj->to_String . "</table>";

	if ($user_can_edit) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-button-icon-left ui-corner-all' onClick="add_experiment_annotation({eid: $eid});"><span class="ui-icon ui-icon-plus"></span>Add Annotation</span>};
	}

	return $html;
}

sub add_experiment_annotation {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;

	my $exp = $coge->resultset('Experiment')->find($eid);

#	my $groups = $coge->resultset('AnnotationTypeGroup');
#	while( my $group = $groups->next ) {
#		foreach my $type ($group->annotation_types) {
#			push @types, { type_name => $type->name, type_id => $type->id };
#		}		
#	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentView.tmpl' );
	$template->param( ADD_EXPERIMENT_ANNOTATION => 1 );
	$template->param( EID            => $eid );

	my %data;
	$data{title} = 'Add Annotation';
	$data{output} = $template->output;

	return encode_json( \%data );
}

sub add_annotation_to_experiment {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;
	my $type_group = $opts{type_group};
	my $type = $opts{type};
	return 0 unless $type;
	my $annotation = $opts{annotation};
	my $link = $opts{link};
	
#	print STDERR "add_annotation_to_list: $lid $type $annotation $link\n";	
	
	my $group_rs;
	if ($type_group) {
		$group_rs = $coge->resultset('AnnotationTypeGroup')->find( { name => $type_group } );
		
		# Create type group if it doesn't already exist
		if (not $group_rs) {
			$group_rs = $coge->resultset('AnnotationTypeGroup')->create( { name => $type_group } );
#			print STDERR "created annotation_type_group " . $group_rs->id . "\n";
		}
	}
	
	my $type_rs;
	$type_rs = $coge->resultset('AnnotationType')->find( 
		{ name => $type, 
		  annotation_type_group_id => ($group_rs ? $group_rs->id : undef) 
		} );
	
	# Create type if it doesn't already exist
	if (not $type_rs) {
		$type_rs = $coge->resultset('AnnotationType')->create( 
			{ name => $type, 
			  annotation_type_group_id => ($group_rs ? $group_rs->id : undef)
			} );
#		print STDERR "created annotation_type " . $type_rs->id . "\n";
	}
	
#	print STDERR "type_rs.id=" . ($type_rs ? $type_rs->id : 'undef') . "\n";
	
	# Create the annotation
	$coge->resultset('ExperimentAnnotation')->create( 
		{ experiment_id => $eid, 
		  annotation => $annotation, 
		  link => $link,
		  annotation_type_id => $type_rs->id
		} );	
	
	return 1;
}

sub remove_experiment_annotation {
	my %opts  = @_;
	my $eid = $opts{eid};
	return "No experiment ID specified" unless $eid;
	my $eaid = $opts{eaid};
	return "No experiment annotation ID specified" unless $eaid;	
	#return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $ea = $coge->resultset('ExperimentAnnotation')->find( { experiment_annotation_id => $eaid } );
	$ea->delete();
	
	return 1;
}

sub generate_html {
	my $html;
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
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
#	my $box_name = "Experiment List: ";
#	my $list_name = $FORM->param('list_name') || $FORM->param('ln');
#	$box_name .= " $list_name" if $list_name;
#	$box_name .= "<span class=link onclick=window.open('$link');>$link</span>";
#	$template->param( BOX_NAME   => $box_name );

	my ($body) = generate_body();
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
		#		foreach my $ug ( $exp->user_groups ) {
		#			$groups .= "<tr>";
		#
		#			$groups .= "<td>" . $ug->name . ( $ug->description ? ': ' . $ug->description : '' ) . "</td>";
		#			$groups .= "<td>" . $ug->role->name . "</td>";
		#			$groups .= "<td>" . join( ", ", map { $_->name } $ug->role->permissions ) . "</td>";
		#
		#			my @users;
		#			foreach my $user ( $ug->users ) {
		#				push @users, $user->name;
		#			}
		#			$groups .= "<td>" . join( ",<br>", @users ) . "</td>";
		#			my $ugid = $ug->id;
		#			$groups .= qq{<td><span class='ui-button ui-corner-all ui-button-icon-left'><span class="ui-icon ui-icon-trash" onclick="} . "remove_group({ugid: $ugid, expid: $expid});" . qq{"</span></span></td>};
		#
		#			$groups .= "</tr>";
		#		}
	}
	return $groups;
}

sub generate_body {
	my $eid = $FORM->param('eid');
	return "Need a valid experiment id\n" unless $eid;
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentView.tmpl' );
	$template->param( MAIN => 1 );
	$template->param( EXPERIMENT_INFO => get_experiment_info(eid=>$eid) );
	$template->param( EXPERIMENT_ANNOTATIONS => get_experiment_annotations(eid => $eid) );
	return $template->output;

#	my $EXPID  = $expid;
#	my $exp;
#	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentView.tmpl' );
#
#	# Create make private/public button
#	my $access_button;
#	my $edit_button;
#	if (1)
#	{    #$USER->is_owner(dsg=>$dsgid) || $USER->is_admin) {
#		if ( $exp->restricted )
#		{
#			$access_button = qq{<span class="ui-button ui-corner-all ui-button-go" onClick="make_experiment_public('$EXPID')">Make Experiment Public</span>};
#		}
#		else
#		{
#			$access_button = qq{<span class="ui-button ui-corner-all ui-button-go" onClick="make_experiment_private('$EXPID')">Make Experiment Private</span>};
#		}
#		$edit_button = qq{<span class="ui-button ui-corner-all ui-button-go" onClick="edit_experiment_info('$EXPID')">Edit Experiment Info</span>};
#	}
#
#	# Build sub-table of types
#	my %types;
#	foreach ( $exp->types )
#	{
#		my $key = $_->name . ( $_->desc ? ': ' . $_->desc : '' );
#		$types{$key}++;
#	}
#	my $type_tbl = '<table class="small"><tbody>';
#	foreach ( sort keys %types )
#	{
#		$type_tbl .= "<tr><td>$_</td></tr>";
#	}
#	$type_tbl .= '</tbody></table>';
#
#	# Build source link
#	my $src  = $exp->source;
#	my $link = $src->link;
#	$link = 'http://' . $link if ( $link and not $link =~ /^http:\/\// );
#	my $source;
#	$source .= "<span class='link' onclick=window.open('$link');>" if ($link);
#	$source .= $src->name . ( $src->desc ? ': ' . $src->desc : '' );
#	$source .= "</span>" if ($link);
#
#	# Build sub-table of annotations
#	%types = ();
#	foreach ( $exp->annotations )
#	{
#		my $type      = $_->annotation_type;
#		my $group     = $type->group;
#		my $groupname = ( $group ? $group->name : '' );
#		my $typename  = $type->name;
#		push @{ $types{$groupname}{$typename} }, $_->annotation;
#	}
#	my $annot_tbl;
#	foreach my $groupname ( sort keys %types )
#	{
#		my $first1 = 1;
#		foreach my $typename ( sort keys %{ $types{$groupname} } )
#		{
#			my $first2 = 1;
#			foreach my $annot ( sort @{ $types{$groupname}{$typename} } )
#			{
#				$annot_tbl .= '<tr>';
#				$annot_tbl .= '<td>' . ( $first1 ? $groupname : '' ) . '</td>'
#					if ($groupname);
#				$annot_tbl .= '<td>' . ( $first2 ? $typename : '' ) . '</td>';
#				$annot_tbl .= '<td></td>' if ( not $groupname );
#				$annot_tbl .= '<td>' . $annot . '</td></tr>';
#				$first1 = $first2 = 0;
#			}
#		}
#	}
#	# Build sub-table of user groups
#	my $groups = get_groups( expid => $EXPID );
#
#	# Populate the template
#	$template->param( EXP_NAME      => $exp->name );
#	$template->param( EXP_DESC      => $exp->desc );
#	$template->param( EXP_TYPE      => $type_tbl );
#	$template->param( EXP_SOURCE    => $source );
#	$template->param( EXP_VER       => $exp->version );
#	$template->param( ACCESS_BUTTON => $access_button );
#	$template->param( EDIT_BUTTON   => $edit_button );
#	$template->param( EXP_ANNOT     => $annot_tbl );
#	$template->param( GROUPS        => $groups );
#
#	return $template->output;
}

sub get_experiment_info {
	my %opts = @_;
	my $eid = $opts{eid};
	my ($exp) = $coge->resultset('Experiment')->find($eid);
	return "Unable to find an entry for $eid" unless $exp;

	my $html = $exp->annotation_pretty_print_html;
	
	if ($USER->is_admin || $USER->is_owner_editor(experiment => $eid)) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="edit_experiment_info({eid: '$eid'});">Edit Experiment Info</span>};
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="add_experiment_type({eid: '$eid'});">Add Experiment Type</span>};
	}
	
	if ($USER->is_admin || $USER->is_owner(experiment => $eid)) {
		if ( $exp->restricted ) {
			$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="make_experiment_public({eid: '$eid'});">Make Experiment Public</span>};
		}
		else {
			$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="make_experiment_private({eid: '$eid'});">Make Experiment Private</span>};
		}
	}	
	
	return $html;
}

sub search_annotation_types {
	my %opts = @_;
	my $type_group = $opts{type_group};
	my $search_term = $opts{search_term};
#	print STDERR "search_annotation_types: $search_term $type_group\n";
	return '' unless $search_term;
	
	$search_term = '%'.$search_term.'%';
	
	my $group;
	if ($type_group) {
		$group = $coge->resultset('AnnotationTypeGroup')->find( { name => $type_group } );
	}
	
	my @types;
	if ($group) {
#		print STDERR "type_group=$type_group " . $group->id . "\n";
		@types = $coge->resultset("AnnotationType")->search(
			\[ 'annotation_type_group_id = ? AND (name LIKE ? OR description LIKE ?)', 
			['annotation_type_group_id', $group->id], ['name', $search_term ], ['description', $search_term] ]);
	}
	else {
		@types = $coge->resultset("AnnotationType")->search(
			\[ 'name LIKE ? OR description LIKE ?', 
			['name', $search_term ], ['description', $search_term] ]);
	}
	
	my %unique;
	map { $unique{$_->name}++ } @types;
	return encode_json( [ sort keys %unique ] );	
}

sub get_annotation_type_groups {
	#my %opts = @_;
	
	my %unique;
	
	my $rs = $coge->resultset('AnnotationTypeGroup');
	while (my $atg = $rs->next) {
		$unique{$atg->name}++;
	}	
	
	return encode_json( [ sort keys %unique ] );
}

