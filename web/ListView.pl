#! /usr/bin/perl -w

use strict;
use CGI;

#use CGI::Ajax;
use JSON::XS;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use Sort::Versions;

#use URI::Escape;
use Data::Dumper;
use File::Path;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;

no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME
  $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION
  $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL
  );
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
$TEMPDIR     = $P->{TEMPDIR} . "ListView/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "ListView/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

%FUNCTION = (
	gen_html          	       => \&gen_html,
	get_list_info     	       => \&get_list_info,
	get_list_contents	       => \&get_list_contents,
	edit_list_info    	       => \&edit_list_info,
	update_list_info  	       => \&update_list_info,
	make_list_public  	       => \&make_list_public,
	make_list_private	       => \&make_list_private,
	add_list_items		       => \&add_list_items,
	add_item_to_list           => \&add_item_to_list,
	remove_list_item	       => \&remove_list_item,
	get_list_annotations       => \&get_list_annotations,
	add_list_annotation        => \&add_list_annotation,
	add_annotation_to_list     => \&add_annotation_to_list,
	remove_list_annotation     => \&remove_list_annotation,
	search_genomes             => \&search_genomes,
	search_experiments         => \&search_experiments,
	search_features            => \&search_features,
	search_lists               => \&search_lists,
	search_annotation_types    => \&search_annotation_types,
	get_annotation_type_groups => \&get_annotation_type_groups,
	delete_list                => \&delete_list,
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
	my $html;
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP       => '/wiki/index.php?title=ListView' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );	
	$template->param( TITLE      => qq{Managing Data} );
	$template->param( PAGE_TITLE => qq{ListView} );
	$template->param( LOGO_PNG   => "ListView-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
	$template->param( ADJUST_BOX => 1 );

	#	$template->param( BOX_NAME	 => $name . " list" );
	$html .= $template->output;
}

sub gen_body {
	my $lid = $FORM->param('lid');
	return "Must have valid list id\n" unless ($lid);
	return "Access denied\n" unless ($USER->has_access(list=>$lid));

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ListView.tmpl' );
	$template->param( PAGE_NAME        => $FORM->url );
	$template->param( MAIN             => 1 );
	$template->param( LIST_INFO        => get_list_info( lid => $lid ) );
	$template->param( LIST_ANNOTATIONS => get_list_annotations( lid => $lid ) );
	$template->param( LIST_CONTENTS    => get_list_contents( lid => $lid ) );
	$template->param( LID              => $lid );
	$template->param( ADMIN_AREA       => 1 ) if $USER->is_admin;
	#$template->param( BOX_NAME        => ... ); #FIXME

	my $open;
	$open = $FORM->param('open') if defined $FORM->param('open');
	my $box_open = $open ? 'true' : 'false';
	$template->param( EDIT_BOX_OPEN => $box_open );

	return $template->output;
}

sub get_list_info {
	my %opts = @_;
	my $lid  = $opts{lid};
	return "Must have valid list id\n" unless ($lid);
	return "Access denied\n" unless ($USER->has_access(list=>$lid));
	
	my ($list) = $coge->resultset('List')->find($lid);
	return "List id$lid does not exist.<br>" .
			"Click <a href='Lists.pl'>here</a> to view a table of all lists." unless ($list);	

	my $html = "List Info:<br>" . $list->annotation_pretty_print_html();

	if ($USER->is_admin || $USER->is_owner_editor(list => $lid)) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="edit_list_info({lid: '$lid'});">Edit List Info</span>};
	}
	
	if ($USER->is_admin || $USER->is_owner(list => $lid)) {
		if ( $list->restricted ) {
			$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="make_list_public({lid: '$lid'});">Make List Public</span>};
		}
		else {
			$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="make_list_private({lid: '$lid'});">Make List Private</span>};
		}
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="dialog_delete_list({lid: '$lid'});">Delete List</span>};
	}

	return $html;
}

sub get_list_types {
	my $current_type_id = shift;
	
	my @types;
	foreach my $type ( $coge->resultset('ListType')->all() ) {
		next if ($type->name =~ /owner/i); # reserve this type for system-created lists
		my $name = $type->name . ($type->description ? ": " . $type->description : '');
		my $selected = '';
		$selected = 'selected="selected"' if ($type->id == $current_type_id);
		push @types, { TID => $type->id, NAME => $name, TYPE_SELECTED => $selected };
	}
	return \@types;
}

sub edit_list_info {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;

	my $list = $coge->resultset('List')->find($lid);
	my $desc = ( $list->description ? $list->description : '' );

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ListView.tmpl' );
	$template->param( EDIT_LIST_INFO => 1 );
	$template->param( LID            => $lid );
	$template->param( NAME           => $list->name );
	$template->param( DESC           => $desc );
	$template->param( TYPE_LOOP      => get_list_types($list->type->id) );

	my %data;
	$data{title} = 'Edit List Info';
	$data{name}   = $list->name;
	$data{desc}   = $desc;
	$data{output} = $template->output;

	return encode_json( \%data );
}

sub update_list_info {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;
	my $name = $opts{name};
	return 0 unless $name;
	my $desc = $opts{desc};
	my $type = $opts{type};

	my $list = $coge->resultset('List')->find($lid);
	$list->name($name);
	$list->description($desc) if $desc;
	$list->list_type_id($type);
	$list->update;

	return 1;
}

sub make_list_public {
	my %opts  = @_;
	my $lid = $opts{lid};
	return "No LID specified" unless $lid;
	#return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $list = $coge->resultset('List')->find($lid);
	$list->restricted(0);
	$list->update;

	return 1;
}

sub make_list_private {
	my %opts  = @_;
	my $lid = $opts{lid};
	return "No LID specified" unless $lid;
	#return "Permission denied." unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	my $list = $coge->resultset('List')->find($lid);
	$list->restricted(1);
	$list->update;

	return 1;
}

sub linkify {
	my ($link, $desc) = @_;
	return "<span class='link' onclick=\"window.open('$link')\">" . $desc . "</span>";
}

sub get_list_annotations {
	my %opts = @_;
	my $lid  = $opts{lid};
	return "Must have valid list id\n" unless ($lid);
	return "Access denied\n" unless ($USER->has_access(list=>$lid));
	
	my ($list) = $coge->resultset('List')->find($lid);
	return unless $list;
	
	my $user_can_edit = !$list->locked && ($USER->is_admin || $USER->is_owner_editor(list => $lid));
	
	my $anno_obj = new CoGe::Accessory::Annotation( Type => 'anno' );
	$anno_obj->Val_delimit("\n");
	$anno_obj->Add_type(0);
	$anno_obj->String_end("\n");
	
	my %groups;
	foreach my $a ( $list->annotations ) {
		my $group = (defined $a->type->group ? $a->type->group->name . ':' . $a->type->name : $a->type->name);
		push @{$groups{$group}}, $a;
	}
	
	foreach my $group (sort keys %groups) {
		foreach my $a ( sort {$a->id <=> $b->id} @{$groups{$group}} ) { #( $list->annotations ) {
			my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . $group . "</span>" );
			$anno_type->Type_delimit(": <td class=\"data5\">");
			$anno_type->Val_delimit("<br>");
			my $a_info = ($a->link ? linkify($a->link, $a->info) : $a->info);
			$a_info .= ($user_can_edit ? "<span onClick=\"remove_list_annotation({lid: '$lid', laid: '" . $a->id . "'});\" class=\"link ui-icon ui-icon-trash\"></span>" : '');
			$anno_type->add_Annot($a_info);
			$anno_obj->add_Annot($anno_type);
		}
	}
	my $html = "<table cellpadding=0 class='ui-widget-content ui-corner-all small'>" . $anno_obj->to_String . "</table>";

	if ($user_can_edit) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-button-icon-left ui-corner-all' onClick="add_list_annotation({lid: $lid});"><span class="ui-icon ui-icon-plus"></span>Add Annotation</span>};
	}

	return 'List Annotations:<br> ' . $html;
}

sub add_list_annotation {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;

	my $list = $coge->resultset('List')->find($lid);

#	my $groups = $coge->resultset('AnnotationTypeGroup');
#	while( my $group = $groups->next ) {
#		foreach my $type ($group->annotation_types) {
#			push @types, { type_name => $type->name, type_id => $type->id };
#		}		
#	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ListView.tmpl' );
	$template->param( ADD_LIST_ANNOTATION => 1 );
	$template->param( LID => $lid );
	$template->param( DEFAULT_TYPE => 'note' );

	my %data;
	$data{title} = 'Add Annotation';
	$data{output} = $template->output;

	return encode_json( \%data );
}

sub add_annotation_to_list {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;
	my $type_group = $opts{type_group};
	my $type = $opts{type};
	return 0 unless $type;
	my $annotation = $opts{annotation};
	my $link = $opts{link};
#	print STDERR "add_annotation_to_list: $lid $type $annotation $link\n";	
	
	$link =~ s/^\s+//;
	$link = 'http://' . $link if (not $link =~ /^(\w+)\:\/\//);

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
	$coge->resultset('ListAnnotation')->create( 
		{ list_id => $lid, 
		  annotation => $annotation, 
		  link => $link,
		  annotation_type_id => $type_rs->id
		} );	
	
	return 1;
}

sub remove_list_annotation {
	my %opts  = @_;
	my $lid = $opts{lid};
	return "No list ID specified" unless $lid;
	my $laid = $opts{laid};
	return "No list annotation ID specified" unless $laid;	
	#return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $list = $coge->resultset('List')->find($lid);
	return 0 if ($list->locked);
	
	my $la = $coge->resultset('ListAnnotation')->find( { list_annotation_id => $laid } );
	$la->delete();
	
	return 1;
}

sub get_list_contents {
	my %opts = @_;
	my $lid  = $opts{lid};
	return "Must have valid list id\n" unless ($lid);
	return "Access denied\n" unless ($USER->has_access(list=>$lid));
	
	my ($list) = $coge->resultset('List')->find($lid);
	return unless $list;
	
	my $user_can_edit = !$list->locked && ($USER->is_admin || $USER->is_owner_editor(list => $lid));
	
	my $anno_obj = new CoGe::Accessory::Annotation( Type => 'anno' );
	$anno_obj->Val_delimit("\n");
	$anno_obj->Add_type(0);
	$anno_obj->String_end("\n");
	
	# FIXME: if objects had a common superclass these three loops could be combined into one
	foreach my $genome ( sort genomecmp $list->genomes ) {
		my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Genomes" . "</span>" );
		$anno_type->Type_delimit(": <td class=\"data5\">");
		$anno_type->Val_delimit("<br>");
		my $gid = $genome->id;
		my $genome_info = $genome->info;
		#FIXME hardcoded value for item_type below
		$genome_info = qq{<span class=link onclick="window.open('OrganismView.pl?dsgid=$gid')">} . $genome_info . ($user_can_edit ? "</span><span onClick=\"remove_list_item({lid: '$lid', item_type: '2', item_id: '$gid'});\" class=\"link ui-icon ui-icon-trash\">" : '') . "</span>";
		$anno_type->add_Annot($genome_info);
		$anno_obj->add_Annot($anno_type);
	}
	foreach my $experiment (sort experimentcmp $list->experiments ) {
		my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Experiments" . "</span>" );
		$anno_type->Type_delimit(": <td class=\"data5\">");
		$anno_type->Val_delimit("<br>");
		my $expid = $experiment->id;
		my $experiment_info = $experiment->info;
		#FIXME hardcoded value for item_type below
		$experiment_info = qq{<span class=link onclick="window.open('ExperimentView.pl?eid=$expid')">} . $experiment_info . ($user_can_edit ? "</span><span onClick=\"remove_list_item({lid: '$lid', item_type: '3', item_id: '$expid'});\" class=\"link ui-icon ui-icon-trash\">" : '') . "</span>";
		$anno_type->add_Annot($experiment_info);
		$anno_obj->add_Annot($anno_type);
	}
	foreach my $feat ( sort featurecmp $list->features ) {
		my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Features" . "</span>" );
		$anno_type->Type_delimit(": <td class=\"data5\">");
		$anno_type->Val_delimit("<br>");
		my $fid = $feat->id;
		my ($feat_info) = $feat->info;
		#FIXME hardcoded value for item_type below
		$feat_info = qq{<span class=link onclick="window.open('FeatView.pl?fid=$fid')">} . $feat_info . ($user_can_edit ? "</span><span onClick=\"remove_list_item({lid: '$lid', item_type: '4', item_id: '$fid'});\" class=\"link ui-icon ui-icon-trash\">" : '') . "</span>";
		$anno_type->add_Annot($feat_info);
		$anno_obj->add_Annot($anno_type);
	}
	foreach my $list ( sort listcmp $list->lists ) {
		my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Lists" . "</span>" );
		$anno_type->Type_delimit(": <td class=\"data5\">");
		$anno_type->Val_delimit("<br>");
		my $child_lid = $list->id;
		my ($list_info) = $list->info;
		#FIXME hardcoded value for item_type below
		$list_info = qq{<span class=link onclick="window.open('ListView.pl?lid=$child_lid')">} . $list_info . ($user_can_edit ? "</span><span onClick=\"remove_list_item({lid: '$lid', item_type: '1', item_id: '$child_lid'});\" class=\"link ui-icon ui-icon-trash\">" : '') . "</span>";
		$anno_type->add_Annot($list_info);
		$anno_obj->add_Annot($anno_type);
	}	
	
	my $html = "<table cellpadding=0 class='ui-widget-content ui-corner-all small'>" . $anno_obj->to_String . "</table>";

	if ($user_can_edit) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-button-icon-left ui-corner-all' onClick="add_list_items({lid: $lid});"><span class="ui-icon ui-icon-plus"></span>Add Items</span>};
	}

	return 'List Contents:<br>' . $html;
}

sub add_list_items {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;

	my $list = $coge->resultset('List')->find($lid);
	my $desc = ( $list->description ? $list->description : '' );

	#
	# Build list of "My Stuff"
	#

	# Experiments
	my @available_items;
	foreach my $exp (sort experimentcmp $USER->experiments) {
		push @available_items, { item_name => 'experiment: ' . $exp->info, 
								 item_spec => 3 . ':' . $exp->id }; #FIXME magic number for item_type
	}
	
	# Genomes
	foreach my $g (sort genomecmp $USER->genomes) {
		push @available_items, { item_name => 'genome: ' . $g->info, 
								 item_spec => 2 . ':' . $g->id }; #FIXME magic number for item_type
	}
	
	# Features
	foreach my $f (sort featurecmp $USER->features) {
		push @available_items, { item_name => 'feature: ' . $f->info, 
								 item_spec => 4 . ':' . $f->id }; #FIXME magic number for item_type
	}
	
	# Lists
	foreach my $l (sort listcmp $USER->lists) {
		next if ($l->id == $lid); # can't add a list to itself!
		next if ($l->locked); # exclude user's master list
		push @available_items, { item_name => 'list: ' . $l->info, 
								 item_spec => 1 . ':' . $l->id }; #FIXME magic number for item_type
	}	
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ListView.tmpl' );
	$template->param( ADD_LIST_ITEMS 	  => 1 );
	$template->param( LID            	  => $lid );
	$template->param( NAME           	  => $list->name );
	$template->param( DESC           	  => $desc );
	$template->param( AVAILABLE_MY_LIST_LOOP => \@available_items );

	# Setup dialog title and data fields
	my %data;
	$data{title}  = 'Add Items to List';
	$data{name}   = $list->name;
	$data{desc}   = $desc;
	$data{items}  = \@available_items;
	$data{output} = $template->output;
	
	return encode_json( \%data );
}

sub add_item_to_list {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;
	my $item_spec = $opts{item_spec};
	return 0 unless $item_spec;
#	print STDERR "$lid $item_spec\n";	
	
	my ($item_type, $item_id) = split(/:/, $item_spec);
	$coge->resultset('ListConnector')->create( { parent_id => $lid, child_id => $item_id, child_type => $item_type } );	
	
	return 1;
}

sub remove_list_item {
	my %opts  = @_;
	my $lid = $opts{lid};
	return "No LID specified" unless $lid;
	#return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $list = $coge->resultset('List')->find($lid);
	return 0 if ($list->locked);
	
	my $item_type = $opts{item_type};
	my $item_id = $opts{item_id};
	
	my $lc = $coge->resultset('ListConnector')->find( { parent_id => $lid, child_id => $item_id, child_type =>$item_type } );
	$lc->delete();
	
	return 1;
}

sub search_genomes {
	my %opts = @_;
	my $lid = $opts{lid};
	my $search_term = $opts{search_term};
#	print STDERR "$lid $search_term\n";
	return 0 unless $lid;	
	return '' unless $search_term and length $search_term > 0;
	
	$search_term = '%'.$search_term.'%';
	
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->genomes;

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
	
	my $html;
	foreach my $g (sort genomecmp @genomes) {
		my $disable = $exists{$g->id} ? "disabled='disabled'" : '';
		my $item_spec = 2 . ':' . $g->id; #FIXME magic number for item_type
		$html .= "<option $disable value='$item_spec'>" . $g->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return $html;
}

sub search_experiments {
	my %opts = @_;
	my $lid = $opts{lid};
	my $search_term = $opts{search_term};
#	print STDERR "$lid $search_term\n";
	return 0 unless $lid;
	return '' unless $search_term;
	
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->experiments;
	
	# Get user's private experiments
	my @experiments;
	foreach ($USER->experiments(restricted => 1)) {
		push @experiments, $_ if ($_->name =~ /$search_term/i or $_->description =~ /$search_term/i);
	}
	
	# Get all public experiments
	$search_term = '%'.$search_term.'%';
	push @experiments, $coge->resultset("Experiment")->search(
		\[ 'restricted=? AND (name LIKE ? OR description LIKE ?)', 
		['restricted', 0], ['name', $search_term ], ['description', $search_term] ]);
	
	my $html;
	foreach my $exp (sort experimentcmp @experiments) {
		my $disable = $exists{$exp->id} ? "disabled='disabled'" : '';
		my $item_spec = 3 . ':' . $exp->id; #FIXME magic number for item_type
		$html .= "<option $disable value='$item_spec'>" . $exp->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return $html;
}

sub search_features {
	my %opts = @_;
	my $lid = $opts{lid};
	my $search_term = $opts{search_term};
#	print STDERR "$lid $search_term\n";
	return 0 unless $lid;
	return '' unless $search_term;
	
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->features;

	# Fulltext search copied from FeatView.pl
    my @fnames;
    push @fnames, $coge->resultset('FeatureName')->search(name => $search_term);
	unless (@fnames) {
		push @fnames, $coge->resultset('FeatureName')->search_literal('MATCH(me.name) AGAINST (?)', $search_term);
	}
	return  "<option disabled='disabled'>" . @fnames . " results, please refine your search.</option>" if (@fnames > 1000);

	my $html;
	my %seen;
	foreach my $f (sort featurecmp @fnames) {
		next if ($seen{$f->feature_id}++);
		my $disable = $exists{$f->feature_id} ? "disabled='disabled'" : '';
		my $item_spec = 4 . ':' . $f->feature_id; #FIXME magic number for item_type
		$html .= "<option $disable value='$item_spec'>" . $f->feature->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return $html;
}

sub search_lists { # list of lists
	my %opts = @_;
	my $lid = $opts{lid};
	my $search_term = $opts{search_term};
#	print STDERR "$lid $search_term\n";
	return 0 unless $lid;
	return '' unless $search_term;
	
	$search_term = '%'.$search_term.'%';
	
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->lists;
	
	my $group_str = join(',', map { $_->id } $USER->groups);
	
	my @lists = $coge->resultset("List")->search_literal(
		"locked=0 AND (restricted=0 OR user_group_id IN ( $group_str )) AND (name LIKE '$search_term' OR description LIKE '$search_term')");
	
	my $html;
	foreach my $l (sort listcmp @lists) {
		next if ($l->id == $lid); # can't add a list to itself!
		my $disable = $exists{$l->id} ? "disabled='disabled'" : '';
		my $item_spec = 1 . ':' . $l->id; #FIXME magic number for item_type
		$html .= "<option $disable value='$item_spec'>" . $l->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
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

sub genomecmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->type->id <=> $b->type->id || $a->name cmp $b->name || $b->id cmp $a->id
}

sub experimentcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->name cmp $b->name || $b->id cmp $a->id
}

sub featurecmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	$a->name cmp $b->name
}

sub listcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	$a->name cmp $b->name
}

sub delete_list {
	my %opts  = @_;
	my $lid = $opts{lid};
	return "No LID specified" unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return "Cannot find list $lid\n" unless $list;
	
	return 0 unless ($USER->is_admin or $USER->is_owner(list => $lid));
	
	if ( $list->locked && !$USER->is_admin ) {
		return "This is a locked list.  Admin permission is needed to modify.";
	}
	
	$list->delete();
	
	return 1;
}
