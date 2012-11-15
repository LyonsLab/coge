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

use vars qw( $P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE 
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

$PAGE_TITLE = "ExperimentView";

$TEMPDIR = $P->{TEMPDIR} . "$PAGE_TITLE/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "$PAGE_TITLE/";

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
	generate_html           => \&generate_html,
	remove_group            => \&remove_group,
	get_groups              => \&get_groups,
	get_experiment_info     => \&get_experiment_info,
	edit_experiment_info    => \&edit_experiment_info,
	update_experiment_info  => \&update_experiment_info,
	get_sources             => \&get_sources,
	make_experiment_public  => \&make_experiment_public,
	make_experiment_private => \&make_experiment_private,
	add_type_to_experiment  => \&add_type_to_experiment,
	get_experiment_types    => \&get_experiment_types,
	get_type_description    => \&get_type_description,
	remove_experiment_type  => \&remove_experiment_type,
	get_annotations			=> \&get_annotations,
	add_annotation			=> \&add_annotation,
	update_annotation		=> \&update_annotation,
	remove_annotation		=> \&remove_annotation,
	get_annotation			=> \&get_annotation,
	search_annotation_types	=> \&search_annotation_types,
	get_annotation_type_groups => \&get_annotation_type_groups,
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
		print $FORM->header, generate_html();
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

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( EDIT_EXPERIMENT_INFO => 1 );
	$template->param( EID            => $eid );
	$template->param( NAME           => $exp->name );
	$template->param( DESC           => $desc );
	$template->param( SOURCE         => $exp->source->name );
	$template->param( SOURCE_ID      => $exp->source->id );
	$template->param( VERSION        => $exp->version );

	my %data;
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
	my $source_id = $opts{source_id};
	my $version = $opts{version};

	my $exp = $coge->resultset('Experiment')->find($eid);
	$exp->name($name);
	$exp->description($desc) if $desc;
	$exp->version($version);
	$exp->data_source_id($source_id);
	$exp->update;

	return 1;
}

sub get_sources {
	#my %opts = @_;
	
	my @sources;
	foreach ($coge->resultset('DataSource')->all()) {
		push @sources, { 'label' => $_->name, 'value' => $_->id };
	}
	
	return encode_json(\@sources);
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

sub add_type_to_experiment {
	my %opts = @_;
	my $eid  = $opts{eid};
	return 0 unless $eid;
	my $name = $opts{name};
	return 0 unless $name;
	my $description = $opts{description};
	
	my $type = $coge->resultset('ExperimentType')->find( { name => $name, description => $description } );
	
	if ($type) {
		# If type exists, check if already assigned to this experiment
		foreach ($type->experiment_type_connectors) {
			return 1 if ($_->experiment_id == $eid);
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
	return "<span class='link' onclick=\"window.open('$link')\">$desc</span>";
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

sub get_annotations {
	my %opts = @_;
	my $eid  = $opts{eid};
	return "Must have valid experiment id\n" unless ($eid);
	my $exp = $coge->resultset('Experiment')->find($eid);
	return "Access denied\n" unless ($USER->has_access(experiment=>$eid) || !$exp->restricted);
	

	return unless $exp;
	
	my $user_can_edit = ($USER->is_admin || $USER->is_owner_editor(experiment => $eid));
	
	my %groups;
	foreach my $a ( $exp->annotations ) {
		my $group = (defined $a->type->group ? $a->type->group->name . ':' . $a->type->name : $a->type->name);
		push @{$groups{$group}}, $a;
	}
	
	my $html = '<table id="experiment_annotation_table" class="ui-widget-content ui-corner-all small" style="max-width:800px;overflow:hidden;word-wrap:break-word;border-spacing:0;border-collapse:collapse;"><thead style="display:none"></thead><tbody>';
	foreach my $group (sort keys %groups) {
		my $first = 1;
		foreach my $a ( sort {$a->id <=> $b->id} @{$groups{$group}} ) {
			$html .= "<tr style='vertical-align:top;'>";
			$html .= "<th align='right' class='title5' style='padding-right:10px;white-space:nowrap;font-weight:normal;background-color:white;' rowspan=" . @{$groups{$group}} . ">$group:</th>" if ($first-- > 0);
			
			#$html .= '<td>';
			my $image_link = ($a->image ? 'image.pl?id=' . $a->image->id : '');
			my $image_info = ($a->image ? "<a href='$image_link' target='_blank' title='click for full-size image'><img height='40' width='40' src='$image_link' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='float:left;padding:1px;border:1px solid lightgray;margin-right:5px;'></a>" : '');
			#$html .= $image_info if $image_info;
			#$html .= "</td>";

			$html .= "<td class='data5'>".$image_info.$a->info.'</td>';
			$html .= '<td style="padding-left:5px;">';
			$html .= linkify($a->link, 'Link') if $a->link;
			$html .= '</td>';
			if ($user_can_edit) {
				my $aid = $a->id;
				$html .= '<td style="padding-left:20px;white-space:nowrap;">' . 
						 "<span onClick=\"edit_annotation_dialog($aid);\" class='link ui-icon ui-icon-gear'></span>" .
						 "<span onClick=\"\$(this).fadeOut(); remove_annotation($aid);\" class='link ui-icon ui-icon-trash'></span>" .
						 '</td>';
			}
			$html .= '</tr>';
		}
	}
	$html .= '</tbody></table>';

	if ($user_can_edit) {
		$html .= qq{<span onClick="add_annotation_dialog();" style="font-size: .75em" class='ui-button ui-button-go ui-button-icon-left ui-corner-all'><span class="ui-icon ui-icon-plus"></span>Add Annotation</span>};
	}

	return $html;
}

sub get_annotation {
	my %opts  = @_;
	my $aid = $opts{aid};
	return unless $aid;
	#TODO check user access here
	
	my $ea = $coge->resultset('ExperimentAnnotation')->find($aid);
	return unless $ea;
	
	my $type = '';
	my $type_group = '';
	if ($ea->type) {
		$type = $ea->type->name;
		$type_group = $ea->type->group->name if ($ea->type->group);
	}
	return encode_json({ annotation => $ea->annotation, link => $ea->link, type => $type, type_group => $type_group });
}

sub add_annotation {
	my %opts = @_;
	my $eid  = $opts{parent_id};
	return 0 unless $eid;
	my $type_group = $opts{type_group};
	my $type = $opts{type};
	return 0 unless $type;
	my $annotation = $opts{annotation};
	my $link = $opts{link};
	my $image_filename = $opts{edit_annotation_image};
	my $fh = $FORM->upload('edit_annotation_image');	
#	print STDERR "add_annotation: $eid $type $annotation $link\n";	
	
	# Create the type and type group if not already present
	my $group_rs;
	if ($type_group) {
		$group_rs = $coge->resultset('AnnotationTypeGroup')->find_or_create( { name => $type_group } );
	}
	my $type_rs = $coge->resultset('AnnotationType')->find_or_create( 
		{ name => $type, 
		  annotation_type_group_id => ($group_rs ? $group_rs->id : undef) 
		} );
	
	# Create the image
	my $image;
	if ($fh) {
		read($fh, my $contents, -s $fh);
		$image = $coge->resultset('Image')->create(
			{	filename => $image_filename,
				image => $contents
			}
		);
		return 0 unless $image;
	}
	
	# Create the annotation
	my $annot = $coge->resultset('ExperimentAnnotation')->create( 
		{ experiment_id => $eid, 
		  annotation => $annotation, 
		  link => $link,
		  annotation_type_id => $type_rs->id,
		  image_id => ($image ? $image->id : undef)
		} );
	return 0 unless $annot;	
	
	return 1;
}

sub update_annotation {
	my %opts  = @_;
	my $aid = $opts{aid};
	return unless $aid;
	my $type_group = $opts{type_group};
	my $type = $opts{type};
	return 0 unless $type;
	my $annotation = $opts{annotation};
	my $link = $opts{link};
	my $image_filename = $opts{edit_annotation_image};
	my $fh = $FORM->upload('edit_annotation_image');	
	#TODO check user access here
	
	my $ea = $coge->resultset('ExperimentAnnotation')->find($aid);
	return unless $ea;

	# Create the type and type group if not already present
	my $group_rs;
	if ($type_group) {
		$group_rs = $coge->resultset('AnnotationTypeGroup')->find_or_create( { name => $type_group } );
	}
	my $type_rs;
	my $type_rs = $coge->resultset('AnnotationType')->find_or_create( 
		{ name => $type, 
		  annotation_type_group_id => ($group_rs ? $group_rs->id : undef) 
		} );
	
	# Create the image
	#TODO if image was changed delete previous image
	my $image;
	if ($fh) {
		read($fh, my $contents, -s $fh);
		$image = $coge->resultset('Image')->create(
			{	filename => $image_filename,
				image => $contents
			}
		);
		return 0 unless $image;
	}
	
	$ea->annotation($annotation);
	$ea->link($link);
	$ea->annotation_type_id($type_rs->id);
	$ea->image_id($image->id) if ($image);
	$ea->update;
	
	return;
}

sub remove_annotation {
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
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( PAGE_TITLE => $PAGE_TITLE );
	$template->param( HELP       => '/wiki/index.php?title=' . $PAGE_TITLE );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= ' ' . $USER->last_name if ( $USER->first_name && $USER->last_name );
	$template->param( USER     => $name );
	$template->param( LOGO_PNG => "$PAGE_TITLE-logo.png" );
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
	$template->param( BODY => $body );

	$template->param( ADJUST_BOX => 1 );
	return $template->output;
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
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( MAIN => 1 );
	$template->param( PAGE_NAME => $PAGE_TITLE . '.pl');
	$template->param( EXPERIMENT_INFO => get_experiment_info(eid=>$eid) );
	$template->param( EXPERIMENT_ANNOTATIONS => get_annotations(eid => $eid) );
	$template->param( EID => $eid );
	$template->param( DEFAULT_TYPE => 'note' );
	
	return $template->output;
}

sub get_experiment_info {
	my %opts = @_;
	my $eid = $opts{eid};
	my ($exp) = $coge->resultset('Experiment')->find($eid);
	return "Access denied\n" unless ($USER->has_access(experiment=>$eid) || !$exp->restricted);

	return "Unable to find an entry for $eid" unless $exp;

	my $allow_edit = $USER->is_admin || $USER->is_owner_editor(experiment => $eid);
	my $html = $exp->annotation_pretty_print_html(allow_delete => $allow_edit);
	
	if ($allow_edit) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="edit_experiment_info();">Edit Info</span>};
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="\$('#experiment_type_edit_box').dialog('open');">Add Type</span>};
	}
	
	if ($USER->is_admin || $USER->is_owner(experiment => $eid)) {
		if ( $exp->restricted ) {
			$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="make_experiment_public();">Make Public</span>};
		}
		else {
			$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="make_experiment_private();">Make Private</span>};
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

