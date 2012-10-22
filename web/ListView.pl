#! /usr/bin/perl -w

use strict;
use CGI;

#use CGI::Ajax;
use JSON::XS;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use Sort::Versions;
use List::Util qw(first);

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

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME
  $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION
  $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL $MAX_SEARCH_RESULTS);
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_NAME = 'ListView.pl';

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

$MAX_SEARCH_RESULTS = 4000;

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( cookie_name => $COOKIE_NAME, ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
$link = CoGe::Accessory::Web::get_tiny_link( db => $coge, user_id => $USER->id, page => $PAGE_NAME, url => $link );


%FUNCTION = (
	gen_html          	       	=> \&gen_html,
	get_list_info     	       	=> \&get_list_info,
	get_list_contents	       	=> \&get_list_contents,
	edit_list_info    	       	=> \&edit_list_info,
	update_list_info  	       	=> \&update_list_info,
	make_list_public  	       	=> \&make_list_public,
	make_list_private	       	=> \&make_list_private,
	add_list_items		       	=> \&add_list_items,
	add_item_to_list           	=> \&add_item_to_list,
	remove_list_item	       	=> \&remove_list_item,
	get_list_annotations       	=> \&get_list_annotations,
	add_list_annotation        	=> \&add_list_annotation,
	add_annotation_to_list     	=> \&add_annotation_to_list,
	remove_list_annotation     	=> \&remove_list_annotation,
	search_genomes             	=> \&search_genomes,
	search_experiments         	=> \&search_experiments,
	search_features            	=> \&search_features,
	search_lists               	=> \&search_lists,
	search_annotation_types    	=> \&search_annotation_types,
	get_annotation_type_groups 	=> \&get_annotation_type_groups,
	delete_list                	=> \&delete_list,
	send_to_blast		   		=> \&send_to_blast,
	send_to_featlist		   		=> \&send_to_featlist,
	send_to_msa					=> \&send_to_msa,
	send_to_gevo				=> \&send_to_gevo,
	send_to_synfind				=> \&send_to_synfind,
	send_to_featmap				=> \&send_to_featmap,
	send_to_codeon				=> \&send_to_codeon,
	send_to_fasta				=> \&send_to_fasta,
	send_to_csv			   		=> \&send_to_csv,
	send_to_xls					=> \&send_to_xls,
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
	#$template->param( BOX_NAME	 => $name . " list" );

	return $template->output;
}

sub gen_body {
	my $lid = $FORM->param('lid');
	return "Must have valid list id\n" unless ($lid);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'ListView.tmpl' );
	$template->param( PAGE_NAME        => $FORM->url );
	$template->param( MAIN             => 1 );
	$template->param( LIST_INFO        => get_list_info( lid => $lid ) );
	$template->param( LIST_ANNOTATIONS => get_list_annotations( lid => $lid ) );
	$template->param( LIST_CONTENTS    => get_list_contents( lid => $lid ) );
	$template->param( LID              => $lid );
	$template->param( ADMIN_AREA       => 1 ) if $USER->is_admin;
	#$template->param( BOX_NAME        => ... ); #FIXME

	return $template->output;
}

sub get_list_info {
	my %opts = @_;
	my $lid  = $opts{lid};
	return unless ($lid);
	return unless ($USER->has_access(list=>$lid));
	
	my ($list) = $coge->resultset('List')->find($lid);
	return "List id$lid does not exist.<br>" .
			"Click <a href='Lists.pl'>here</a> to view a table of all lists." unless ($list);	

	my $html = $list->annotation_pretty_print_html();

	if ($USER->is_admin || (not $list->locked && $USER->is_owner_editor(list => $lid))) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-corner-all' onClick="edit_list_info({lid: '$lid'});">Edit List Info</span>};
	}
	
	if ($USER->is_admin || (not $list->locked && $USER->is_owner(list => $lid))) {
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
	return unless $lid;
	#return unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $list = $coge->resultset('List')->find($lid);
	$list->restricted(0);
	$list->update;

	return 1;
}

sub make_list_private {
	my %opts  = @_;
	my $lid = $opts{lid};
	return unless $lid;
	#return unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
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
	return unless ($lid);
	return unless ($USER->has_access(list=>$lid));
	
	my ($list) = $coge->resultset('List')->find($lid);
	return unless $list;
	
	my $user_can_edit = $USER->is_admin || (!$list->locked && $USER->is_owner_editor(list => $lid));
	
	my $anno_obj = new CoGe::Accessory::Annotation( Type => 'anno' );
	$anno_obj->Val_delimit("\n");
	$anno_obj->Add_type(0);
	$anno_obj->String_end("\n");
	
	my %groups;
	foreach my $a ( $list->annotations ) {
		my $group = (defined $a->type->group ? $a->type->group->name . ':' . $a->type->name : $a->type->name);
		push @{$groups{$group}}, $a;
	}
	
	my $html = '<table id="list_annotation_table" cellpadding=0 class="ui-widget-content ui-corner-all small" style="max-width:400px;overflow:hidden;word-wrap:break-word;"><thead style="display:none"></thead><tbody>';
	foreach my $group (sort keys %groups) {
		my $first = 1;
		foreach my $a ( sort {$a->id <=> $b->id} @{$groups{$group}} ) {
			$html .= "<tr valign='top'>" . ($first-- > 0 ? "<th align='right' class='title5' nowrap='true' rowspan=" . @{$groups{$group}} . " style='font-weight:normal;background-color:white'>$group:</th>" : '');
			
			my $image_link = ($a->image ? 'image.pl?id=' . $a->image->id : '');
			my $a_info = ($a->image ? "<a href='$image_link' target='_blank'><img height=20 width=20 src='$image_link' style='vertical-align:text-top;'></a>" : '');
			$a_info .= ' ' . ($a->link ? linkify($a->link, $a->info) : $a->info);
			$html .= "<td class='data5'>$a_info</td>";
			
			if ($user_can_edit) {
				$html .= "<td><span onClick=\"\$(this).fadeOut(); remove_list_annotation({lid: '$lid', laid: '" . $a->id . "'});\" class='link ui-icon ui-icon-trash'></span></td>";
			}
			$html .= '</tr>';
		}
	}
	$html .= '</tbody></table>';
	
	if ($user_can_edit) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-button-icon-left ui-corner-all' onClick="add_list_annotation({lid: $lid});"><span class="ui-icon ui-icon-plus"></span>Add Annotation</span>};
	}	

	return $html;
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
	my $image_filename = $opts{edit_annotation_image};
	my $fh = $FORM->upload('edit_annotation_image');
	#return "Image file is too large (>10MB)" if (-s $fh > 10*1024*1024); # FIXME
	
	if ($link) {
		$link =~ s/^\s+//;
		$link = 'http://' . $link if (not $link =~ /^(\w+)\:\/\//);
	}

	my $group_rs;
	if ($type_group) {
		$group_rs = $coge->resultset('AnnotationTypeGroup')->find( { name => $type_group } );
		
		# Create type group if it doesn't already exist
		if (not $group_rs) {
			$group_rs = $coge->resultset('AnnotationTypeGroup')->create( { name => $type_group } );
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
	}
	
	# Create the annotation
	my $la = $coge->resultset('ListAnnotation')->create( 
		{ list_id => $lid, 
		  annotation => $annotation, 
		  link => $link,
		  annotation_type_id => $type_rs->id
		} );

	# Create the image
	if ($fh) {
		read($fh, my $image, -s $fh);
		$coge->resultset('Image')->create(
			{	list_annotation_id => $la->id,
				filename => $image_filename,
				image => $image
			}
		);
	}
	
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
	return 0 if ($list->locked and not $USER->is_admin);
	
	my $la = $coge->resultset('ListAnnotation')->find( { list_annotation_id => $laid } );
	$la->delete();
	
	return 1;
}

sub get_list_contents {
	my %opts = @_;
	my $lid  = $opts{lid};
	return "Must have valid list id\n" unless ($lid);
	
	my $list = $coge->resultset('List')->find($lid);
	return "List id$lid does not exist.<br>" .
			"Click <a href='Lists.pl'>here</a> to view all lists." unless $list;
			
	return "Access denied\n" unless ($USER->has_access(list=>$lid));
		
	my $user_can_edit = $USER->is_admin || (!$list->locked && $USER->is_owner_editor(list => $lid));
	
	my $html;
	my $first = 1;
	$html = '<table id="list_contents_table" class="small ui-widget-content ui-corner-all"><thead style="display:none"></thead><tbody>';

	my $list_types = CoGeX::list_child_types();
	my $genome_count = $list->genomes(count=>1); #EL: moved outside of loop; massive speed improvement due to cost of this call
	foreach my $genome ( sort genomecmp $list->genomes ) {
		
		$html .= "<tr valign='top'>" . ($first-- > 0 ? "<th align='right' class='title5' nowrap='true' rowspan=$genome_count style='font-weight:normal;background-color:white'>Genomes ($genome_count):</th>" : '');
		my $gid = $genome->id;
		$html .= qq{<td class='data5'><span id='genome$gid' class='link' onclick="window.open('OrganismView.pl?dsgid=$gid')">} . $genome->info . "</span></td>";
		if ($user_can_edit) {
			$html .= "<td><span onClick=\"remove_list_item(this, {lid: '$lid', item_type: '".$list_types->{genome}."', item_id: '$gid'});\" class='link ui-icon ui-icon-circle-minus'></span></td>";
		}
		$html .= '</tr>';
	}
	$first = 1;
	my $exp_count = $list->experiments(count=>1); #EL: moved outside of loop; massive speed improvement
	foreach my $experiment (sort experimentcmp $list->experiments ) {

		$html .= "<tr valign='top'>" . ($first-- > 0 ? "<th align='right' class='title5' nowrap='true' rowspan=$exp_count style='font-weight:normal;background-color:white'>Experiments ($exp_count):</th>" : '');
		my $eid = $experiment->id;
		$html .= qq{<td class='data5'><span id='experiment$eid' class='link' onclick="window.open('ExperimentView.pl?eid=$eid')">} . $experiment->info . "</span></td>";
		if ($user_can_edit) {
			$html .= "<td><span onClick=\"remove_list_item(this, {lid: '$lid', item_type: '".$list_types->{experiment}."', item_id: '$eid'});\" class='link ui-icon ui-icon-circle-minus'></span></td>";
		}
		$html .= '</tr>';
	}
	$first = 1;
	my $feat_count = $list->features(count=>1); #EL: moved outside of loop; massive speed improvement
	foreach my $feature (sort featurecmp $list->features ) {

		$html .= "<tr valign='top'>" . ($first-- > 0 ? "<th align='right' class='title5' nowrap='true' rowspan=$feat_count style='font-weight:normal;background-color:white'>Features ($feat_count):</th>" : '');
		my $fid = $feature->id;
		$html .= qq{<td class='data5'><span id='feature$fid' class='link' onclick="window.open('FeatView.pl?fid=$fid')">} . $feature->info . "</span></td>";
		if ($user_can_edit) {
			$html .= "<td><span onClick=\"remove_list_item(this, {lid: '$lid', item_type: '".$list_types->{feature}."', item_id: '$fid'});\" class='link ui-icon ui-icon-circle-minus'></span></td>";
		}
		$html .= '</tr>';
	}
	$first = 1;
	my $list_count = $list->lists(count=>1);
	foreach my $list (sort listcmp $list->lists ) {
		$html .= "<tr valign='top'>" . ($first-- > 0 ? "<th align='right' class='title5' nowrap='true' rowspan=$list_count style='font-weight:normal;background-color:white'>Lists ($list_count):</th>" : '');
		my $child_id = $list->id;
		$html .= qq{<td class='data5'><span id='list$child_id' class='link' onclick="window.open('ListView.pl?lid=$child_id')">} . $list->info . "</span></td>";
		if ($user_can_edit) {
			$html .= "<td><span onClick=\"remove_list_item(this, {lid: '$lid', item_type: '".$list_types->{list}."', item_id: '$child_id'});\" class='link ui-icon ui-icon-circle-minus'></span></td>";
		}
		$html .= '</tr>';
	}		

	$html .= '</tbody></table>';
	
	if ($user_can_edit) {
		$html .= qq{<span style="font-size: .75em" class='ui-button ui-button-go ui-button-icon-left ui-corner-all' onClick="add_list_items({lid: $lid});"><span class="ui-icon ui-icon-plus"></span>Add Items</span>};
		
	}
	
	return $html;
}

sub add_list_items {
	my %opts = @_;
	my $lid  = $opts{lid};
	return 0 unless $lid;

	my $list = $coge->resultset('List')->find($lid);
	my $desc = ( $list->description ? $list->description : '' );

	my $child_types = CoGeX::list_child_types();

	# Experiments

	my @available_items;
	my %exists = map { $_->id => 1 } $list->experiments;
	foreach my $e (sort experimentcmp $USER->experiments) {
		push @available_items, { item_name => 'experiment: ' . $e->info, 
								 item_spec => $child_types->{experiment} . ':' . $e->id,
								 item_disable => ($exists{$e->id} ? "disabled='disabled'" : '') };
	}

	# Genomes
	%exists = map { $_->id => 1 } $list->genomes;
	foreach my $g (sort genomecmp $USER->genomes) {
		push @available_items, { item_name => 'genome: ' . $g->info, 
								 item_spec =>  $child_types->{genome} . ':' . $g->id,
								 item_disable => ($exists{$g->id} ? "disabled='disabled'" : '') };
	}
	# Features
	%exists = map { $_->id => 1 } $list->features;
	foreach my $f (sort featurecmp $USER->features) {
		push @available_items, { item_name => 'feature: ' . $f->info, 
								 item_spec =>  $child_types->{feature} . ':' . $f->id,
								 item_disable => ($exists{$f->id} ? "disabled='disabled'" : '') };
	}
	# Lists
	%exists = map { $_->id => 1 } $list->lists;
	foreach my $l (sort listcmp $USER->lists) {
		next if ($l->id == $lid); # can't add a list to itself!
		next if ($l->locked); # exclude user's master list
		push @available_items, { item_name => 'list: ' . $l->info, 
								 item_spec =>  $child_types->{list} . ':' . $l->id, 
								 item_disable => ($exists{$l->id} ? "disabled='disabled'" : '') };
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
	
	my ($item_type, $item_id) = split(/:/, $item_spec);
	
	my $lc = $coge->resultset('ListConnector')->create( { parent_id => $lid, child_id => $item_id, child_type => $item_type } );	

	return 1;
}

sub remove_list_item {
	my %opts  = @_;
	my $lid = $opts{lid};
	return "No LID specified" unless $lid;
	#return "Permission denied" unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
	
	my $list = $coge->resultset('List')->find($lid);
	return 0 if ($list->locked and not $USER->is_admin);
	
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
	my $timestamp = $opts{timestamp};
#	print STDERR "$lid $search_term $timestamp\n";
	return 0 unless $lid;	

	# Get genomes already in list
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->genomes;

	my %unique;
	my $num_results;
	
	# Try to get all items if blank search term
	if (not $search_term) {
		# Get all genomes
		$num_results = $coge->resultset("Genome")->count;
		if ($num_results < $MAX_SEARCH_RESULTS) {
			my @genomes = $coge->resultset("Genome")->all;
			map { $unique{$_->id} = $_ if (not $_->restricted or $USER->has_access_to_genome($_)) } @genomes;
		}
	}
	# Perform search
	else {
		$search_term = '%'.$search_term.'%';

		# Get all matching organisms
		my @organisms = $coge->resultset("Organism")->search(
			\[ 'name LIKE ? OR description LIKE ?', 
			['name', $search_term ], ['description', $search_term] ]);

		# Get all matching genomes
		my @genomes = $coge->resultset("Genome")->search(
			\[ 'name LIKE ? OR description LIKE ?', 
			['name', $search_term], ['description', $search_term] ]);

		# Combine matching genomes with matching organism genomes, preventing duplicates
		map { $unique{$_->id} = $_ if (not $_->restricted or $USER->has_access_to_genome($_)) } @genomes;
		foreach my $organism (@organisms) {
			map { $unique{$_->id} = $_ if (not $_->restricted or $USER->has_access_to_genome($_)) } $organism->genomes;
		}
		
		$num_results = keys %unique;
	}
	
	# Limit number of results displayed
	if ($num_results > $MAX_SEARCH_RESULTS) {
		return encode_json({
					timestamp => $timestamp, 
					html => "<option disabled='disabled'>$num_results results, please refine your search.</option>"
		});
	}
	
	# Build select options out of results
	my $child_types = CoGeX::list_child_types;
	my $html;
	foreach my $g (sort genomecmp values %unique) {
		my $disable = $exists{$g->id} ? "disabled='disabled'" : '';
		my $item_spec = $child_types->{genome} . ':' . $g->id; 
		$html .= "<option $disable value='$item_spec'>" . $g->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return encode_json({timestamp => $timestamp, html => $html});
}

sub search_experiments {
	my %opts = @_;
	my $lid = $opts{lid};
	my $search_term = $opts{search_term};
	my $timestamp = $opts{timestamp};
#	print STDERR "$lid $search_term\n";
	return 0 unless $lid;
	
	# Get experiments already in list
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->experiments;
	
	my @experiments;
	my $num_results;
	
	# Try to get all items if blank search term
	if (not $search_term) {
		# Get all experiments
		$num_results = $coge->resultset("Experiment")->count;
		if ($num_results < $MAX_SEARCH_RESULTS) {
			@experiments = $coge->resultset("Experiment")->all;
		}
	}
	# Perform search
	else {
		# Get user's private experiments
		foreach ($USER->experiments(restricted => 1)) {
			push @experiments, $_ if ($_->name =~ /$search_term/i or $_->description =~ /$search_term/i);
		}
		
		# Get all public experiments
		$search_term = '%'.$search_term.'%';
		push @experiments, $coge->resultset("Experiment")->search(
			\[ 'restricted=? AND (name LIKE ? OR description LIKE ?)', 
			['restricted', 0], ['name', $search_term ], ['description', $search_term] ]);
			
		$num_results = @experiments;
	}
	
	# Limit number of results displayed
	if ($num_results > $MAX_SEARCH_RESULTS) {
		return encode_json({
					timestamp => $timestamp, 
					html => "<option disabled='disabled'>$num_results results, please refine your search.</option>"
		}); 
	}
	
	# Build select items out of results
	my $child_types = CoGeX::list_child_types;
	my $html;
	foreach my $exp (sort experimentcmp @experiments) {
		my $disable = $exists{$exp->id} ? "disabled='disabled'" : '';
		my $item_spec = $child_types->{experiment} . ':' . $exp->id; 
		$html .= "<option $disable value='$item_spec'>" . $exp->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return encode_json({timestamp => $timestamp, html => $html});
}

sub search_features {
	my %opts = @_;
	my $lid = $opts{lid};
	my $search_term = $opts{search_term};
	my $timestamp = $opts{timestamp};
#	print STDERR "$lid $search_term\n";
	return 0 unless $lid;
	
	# Get lists already in list
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->features;

	my @fnames;
	my $num_results;
	
	# Try to display all items if blank search term
	if (not $search_term) {
		# Get all features
		$num_results = $coge->resultset("FeatureName")->count;
		if ($num_results < $MAX_SEARCH_RESULTS) {
			@fnames = $coge->resultset("FeatureName")->all;
		}
	}
	# Perform search
	else {
		# Fulltext search copied from FeatView.pl
	    push @fnames, $coge->resultset('FeatureName')->search(name => $search_term);
		unless (@fnames) {
			push @fnames, $coge->resultset('FeatureName')->search_literal('MATCH(me.name) AGAINST (?)', $search_term);
		}
		$num_results = @fnames;
	}
	
	# Limit the number of results displayed
	if ($num_results > $MAX_SEARCH_RESULTS) {
		return encode_json({
					timestamp => $timestamp, 
					html => "<option disabled='disabled'>$num_results results, please refine your search.</option>"
		});
	}

	# Build select items out of results
	my $child_types = CoGeX::list_child_types;
	my $html;
	my %seen;
	foreach my $f (sort featurecmp @fnames) {
		next if ($seen{$f->feature_id}++);
		my $disable = $exists{$f->feature_id} ? "disabled='disabled'" : '';
		my $item_spec = $child_types->{feature} . ':' . $f->feature_id; 
		$html .= "<option $disable value='$item_spec'>" . $f->feature->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return encode_json({timestamp => $timestamp, html => $html});
}

sub search_lists { # list of lists
	my %opts = @_;
	my $lid 		= $opts{lid};
	my $search_term	= $opts{search_term};
	my $timestamp	= $opts{timestamp};
#	print STDERR "$lid $search_term $timestamp\n";
	return 0 unless $lid;
	
	# Get lists already in this list
	my $list = $coge->resultset('List')->find($lid);
	my %exists;
	map { $exists{$_->id}++ } $list->lists;

	my @lists;
	my $num_results;
	my $group_str = join(',', map { $_->id } $USER->groups);
	
	# Try to get all items if blank search term
	if (not $search_term) {
		# Get all lists
		my $sql = "locked=0 AND (restricted=0 OR user_group_id IN ( $group_str ))";
		$num_results = $coge->resultset("List")->count_literal($sql);
		if ($num_results < $MAX_SEARCH_RESULTS) {
			@lists = $coge->resultset("List")->search_literal($sql);
		}
	}
	# Perform search
	else {
		# Get public lists and user's private lists	
		$search_term = '%'.$search_term.'%';
		@lists = $coge->resultset("List")->search_literal(
			"locked=0 AND (restricted=0 OR user_group_id IN ( $group_str )) \
			 AND (name LIKE '$search_term' OR description LIKE '$search_term')");
		$num_results = @lists;
	}
	
	# Limit number of results display
	if ($num_results > $MAX_SEARCH_RESULTS) {
		return encode_json({
					timestamp => $timestamp,
					html => "<option disabled='disabled'>$num_results results, please refine your search.</option>"
		});
	}
	
	# Build select items out of results
	my $child_types = CoGeX::list_child_types;
	my $html;
	foreach my $l (sort listcmp @lists) {
		next if ($l->id == $lid); # can't add a list to itself!
		my $disable = $exists{$l->id} ? "disabled='disabled'" : '';
		my $item_spec = $child_types->{list} . ':' . $l->id;
		$html .= "<option $disable value='$item_spec'>" . $l->info . "</option><br>\n";	
	}
	$html = "<option disabled='disabled'>No matching items</option>" unless $html;
	
	return encode_json({timestamp => $timestamp, html => $html});
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

# FIXME these comparison routines are duplicated elsewhere
sub genomecmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	$a->organism->name cmp $b->organism->name || versioncmp($b->version, $a->version) || $a->type->id <=> $b->type->id || $a->name cmp $b->name || $b->id cmp $a->id
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

sub send_to_blast
{
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	#my $list = $coge->resultset('List')->find($lid);
	#return unless $list;		
	
	#my $accn_list = join(',', map { $_->id } $list->genomes);
	#my $url = "CoGeBlast.pl?dsgid=$accn_list";
	my $url = "CoGeBlast.pl?lid=$lid";
	return encode_json({url => $url});
}

sub send_to_featlist
{
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	my $url = "FeatList.pl?lid=$lid";
	return encode_json({url => $url});
}

sub send_to_msa {
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;		
	
	my $accn_list = join('&', map { 'fid=' . $_->id } $list->features);
	my $url = "CoGeAlign.pl?$accn_list";
	return encode_json({url => $url});
}

sub send_to_gevo
{
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;	
	
	my $count = 1;
	my $accn_list = join('&', map { 'fid' . $count++ . '=' . $_->id } $list->features);
		
	my $url = "GEvo.pl?$accn_list";
	$count--;
	return encode_json({alert => "You have exceeded the number of features you can send to GEvo (20 max)."}) if $count > 20;
	$url .= "&num_seqs=$count";
	return encode_json({url => $url});
}

sub send_to_synfind
{
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;	
	
	my $accn_list = join(',', map { $_->id } $list->genomes);
	my $url = "SynFind.pl?dsgid=$accn_list";
	my $first = first { $_->id } $list->features;
	$url .= ";fid=" . $first->id if ($first);
	return encode_json({url => $url});
}

sub send_to_featmap {
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;	
	
	my $accn_list = join('&', map { 'fid=' . $_->id } $list->features);
	my $url = "FeatMap.pl?$accn_list";
	return encode_json({url => $url});
}

sub send_to_codeon {
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;	
	
	my $accn_list = join('::', map { $_->id } $list->features);
	my $url = "CodeOn.pl?fid=$accn_list";
	return encode_json({url => $url});	
}

sub send_to_fasta {
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;

	$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
	my $basename = $cogeweb->basefilename;
	my $file     = $TEMPDIR . "$basename.faa";
	open( OUT, ">$file" );
	foreach my $g ( sort genomecmp $list->genomes ) {
		print OUT $g->fasta;
	}
	close OUT;
	$file =~ s/$TEMPDIR/$TEMPURL/;
	return $file;
}

sub send_to_csv 
{
	my %opts = @_;
	my $lid = $opts{lid};
	return unless $lid;
	
	my $list = $coge->resultset('List')->find($lid);
	return unless $list;	

	$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
	my $basename = $cogeweb->basefilename;
	my $file = "$TEMPDIR/$basename.csv";
#	print STDERR $P->{SERVER} . " $TEMPURL $file\n";
	open( OUT, ">$file" );

	print OUT join( "\t", "CoGe Genome ID", "Name", "Description", "Source", "Provenance", "Sequence Type", "Chr Count", "Length (bp)", "Percent GC", "Percent AT", "Percent N|X", "OrganismView Link" ), "\n";
	foreach my $g ( sort genomecmp $list->genomes ) {
		my $name = $g->name ? $g->name : $g->organism->name;
		my $desc = $g->description ? $g->description : $g->organism->description;
		my ($ds_source) = $g->source;
		my $source = $ds_source->name;
		my $provenance = join( "||", map { $_->name } $g->datasets );
		my $chr_count  = $g->chromosome_count;
		my $length     = $g->length;
		my $type       = $g->type->name;
		my ( $gc, $at, $n ) = $g->percent_gc();
		$at *= 100;
		$gc *= 100;
		print OUT join( "\t", $g->id, $name, $desc, $source, $provenance, $type, $chr_count, $length, $gc, $at, $n, $P->{SERVER} . 'OrganismView.pl?dsgid=' . $g->id ), "\n";
	}
	
	
	close OUT;
	$file =~ s/$TEMPDIR/$TEMPURL/;
	return encode_json({url => $file});	
}

sub send_to_xls {
	my %args      = @_;
	my $accn_list = $args{accn};
	$accn_list =~ s/^,//;
	$accn_list =~ s/,$//;
	$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
	my $basename = $cogeweb->basefilename;
	my $file     = "$TEMPDIR/Excel_$basename.xls";
	my $workbook = Spreadsheet::WriteExcel->new($file);
	$workbook->set_tempdir("$TEMPDIR");
	my $worksheet = $workbook->add_worksheet();
	my $i         = 1;

	$worksheet->write( 0, 0,  "Name" );
	$worksheet->write( 0, 1,  "Description" );
	$worksheet->write( 0, 2,  "Source" );
	$worksheet->write( 0, 3,  "Provenance" );
	$worksheet->write( 0, 4,  "Sequence Type" );
	$worksheet->write( 0, 5,  "Chr Count" );
	$worksheet->write( 0, 6,  "Length (bp)" );
	$worksheet->write( 0, 7,  "Percent GC" );
	$worksheet->write( 0, 8,  "Percent AT" );
	$worksheet->write( 0, 9,  "Percent N|X" );
	$worksheet->write( 0, 10, "OrganismView Link" );

	foreach my $dsgid ( split(/,/, $accn_list) ) {
		my ($dsg) = $coge->resultset("Genome")->find($dsgid);
		next unless $dsg;
		
		my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
		my $desc = $dsg->description ? $dsg->description : $dsg->organism->description;
		my ($ds_source) = $dsg->source;
		my $source = $ds_source->name;
		my $provenance = join( " ", map { $_->name } $dsg->datasets );
		my $length     = $dsg->length;
		my $chr_count  = $dsg->chromosome_count;
		my $type       = $dsg->type->name;
		my ( $gc, $at, $n ) = $dsg->percent_gc();
		$at *= 100;
		$gc *= 100;

		#my ($wgc, $wat) = $dsg->wobble_content;
		#$wat*=100;
		#$wgc*=100;

		$worksheet->write( $i, 0, $name );
		$worksheet->write( $i, 1, $desc );
		$worksheet->write( $i, 2, $source );
		$worksheet->write( $i, 3, $provenance );
		$worksheet->write( $i, 4, $type );
		$worksheet->write( $i, 5, $chr_count );
		$worksheet->write( $i, 6, $length );
		$worksheet->write( $i, 7, $gc . '%' );
		$worksheet->write( $i, 8, $at . '%' );
		$worksheet->write( $i, 9, $n . '%' );
		$worksheet->write( $i, 10, $P->{SERVER} . 'OrganismView.pl?dsgid=' . $dsgid );

		$i++;
	}
	$workbook->close() or die "Error closing file: $!";
	$file =~ s/$TEMPDIR/$TEMPURL/;
	return $file;
}
