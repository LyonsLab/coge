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
use URI::Escape::JavaScript qw(escape unescape);
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use DBIxProfiler;
use File::Path;
use Sort::Versions;
use LWP::UserAgent;
use LWP::Simple;
use HTTP::Status qw(:constants);
use File::Listing;
use File::Copy;
use XML::Simple;
no warnings 'redefine';

use vars qw(
	$P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE 
	$TEMPDIR $BINDIR $USER $DATE $COGEDIR $coge $FORM $URL $TEMPURL $COOKIE_NAME 
	%FUNCTION $MAX_SEARCH_RESULTS
);

$P         = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_TITLE = 'GenomeInfo';

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

$TEMPDIR = $P->{TEMPDIR} . $PAGE_TITLE . '/' . $USER->name . '/';
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;

$MAX_SEARCH_RESULTS = 100;

#$SIG{'__WARN__'} = sub { };    # silence warnings

%FUNCTION = (
	generate_html			=> \&generate_html,
	get_genome_info			=> \&get_genome_info,
	get_genome_data			=> \&get_genome_data,
	get_notebooks			=> \&get_notebooks,
	get_users_with_access	=> \&get_users_with_access,
	get_groups_with_access 	=> \&get_groups_with_access,
	get_available_users		=> \&get_available_users,
	get_available_groups	=> \&get_available_groups,
	edit_genome_info		=> \&edit_genome_info,
	update_genome_info		=> \&update_genome_info,
	update_owner			=> \&update_owner,
	add_user				=> \&add_user,
	remove_user				=> \&remove_user,
	add_to_notebook			=> \&add_to_notebook,
	add_genome_to_notebook	=> \&add_genome_to_notebook,
	remove_genome_from_notebook	=> \&remove_genome_from_notebook,
	get_notebook_preview	=> \&get_notebook_preview,
	search_organisms		=> \&search_organisms,
	search_users			=> \&search_users,
	search_notebooks		=> \&search_notebooks,
);

if ( $FORM->param('jquery_ajax') ) {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {
		die if (not defined $FUNCTION{$fname});
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else {
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
}
else {
	print $FORM->header, "\n", generate_html();
}

sub get_genome_info {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param(
		DO_GENOME_INFO => 1,
		ORGANISM => $genome->organism->name,
		VERSION => $genome->version,
		TYPE => $genome->type->name,
		SOURCE => join(',', map { $_->name } $genome->source),
		RESTRICTED => ($genome->restricted ? 'Yes' : 'No'),
		NAME => $genome->name,
		DESCRIPTION => $genome->description,
	);

	$template->param( GID => $genome->id );

	return $template->output;
}

sub edit_genome_info {
	my %opts = @_;
	my $gid  = $opts{gid};
	return unless ($gid);

	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( 
		EDIT_GENOME_INFO => 1, 
		ORGANISM => $genome->organism->name,
		VERSION => $genome->version,
		TYPE => $genome->type->name,
		SOURCE => join(',', map { $_->name } $genome->source),
		RESTRICTED => $genome->restricted,
		NAME => $genome->name,
		DESCRIPTION => $genome->description 
	);

	$template->param(
		TYPES => get_sequence_types($genome->type->id),
		SOURCES => get_sources()
	);

	return $template->output;
}

sub update_genome_info {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $name = $opts{name};
	my $description = $opts{description};
	my $version = $opts{version};
	my $type_id = $opts{type_id};
	my $restricted = $opts{restricted};
	my $org_name = $opts{org_name};
	my $source_name = $opts{source_name};
	my $timestamp = $opts{timestamp};
	# print STDERR "gid=$gid organism=$org_name version=$version source=$source_name\n";
	return "Error: missing params." unless ($gid and $org_name and $version and $source_name);

	my $genome = $coge->resultset('Genome')->find($gid);
	return "Error: can't find genome." unless ($genome);

	my $organism = $coge->resultset('Organism')->find({name => $org_name});
	return "Error: can't find organism." unless ($organism);

	my $source = $coge->resultset('DataSource')->find({name => $source_name});
	return "Error: can't find source." unless ($source);

	$genome->organism_id($organism->id);
	$genome->name($name);
	$genome->description($description);
	$genome->version($version);
	$genome->genomic_sequence_type_id($type_id);
	$genome->restricted($restricted eq 'true');

	foreach my $ds ($genome->datasets) {
		$ds->data_source_id($source->id);
		$ds->update;
	}

	$genome->update;

	return;
}

sub update_owner {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $user_name = $opts{user_name};
	return unless ($gid and $user_name);

	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my $user = $coge->resultset('User')->find({user_name => $user_name});
	return unless ($user);

	# Remove from current owner list
	if ($genome->owner_list) {
		my $conn = $coge->resultset('ListConnector')->find({parent_id => $genome->owner_list->id, child_id => $gid});
		$conn->delete;
	}

	# Add to user's owner list
	my $conn = $coge->resultset('ListConnector')->create({parent_id => $user->owner_list->id, child_id => $gid});
	return unless ($conn);

	return 1;
}

sub search_organisms {
	my %opts = @_;
	my $search_term = $opts{search_term};
	my $timestamp = $opts{timestamp};
#	print STDERR "$search_term $timestamp\n";
	return unless $search_term;

	# Perform search
	$search_term = '%'.$search_term.'%';
	my @organisms = $coge->resultset("Organism")->search(
		\[ 'name LIKE ? OR description LIKE ?', 
		['name', $search_term ], ['description', $search_term] ]);

	# Limit number of results displayed
	if (@organisms > $MAX_SEARCH_RESULTS) {
		return encode_json({timestamp => $timestamp, items => undef});
	}
	
	my %unique = map { $_->name => 1 } @organisms;
	return encode_json({timestamp => $timestamp, items => [sort keys %unique]});
}

sub search_users {
	my %opts = @_;
	my $search_term = $opts{search_term};
	my $timestamp = $opts{timestamp};
	#print STDERR "$search_term $timestamp\n";
	return unless $search_term;

	# Perform search
	$search_term = '%'.$search_term.'%';
	my @users = $coge->resultset("User")->search(
		\[ 'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?', 
		['user_name', $search_term], ['first_name', $search_term], ['last_name', $search_term] ]);

	# Limit number of results displayed
	# if (@users > $MAX_SEARCH_RESULTS) {
	# 	return encode_json({timestamp => $timestamp, items => undef});
	# }
	
	return encode_json({timestamp => $timestamp, items => [sort map { $_->user_name } @users]});
}

sub get_genome_data {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( 
		DO_GENOME_DATA => 1,
		CHROMOSOME_COUNT => commify($genome->chromosome_count()),
		LENGTH => commify($genome->length),
		GID => $genome->id );
	
	return $template->output;
}

sub get_genome_owner {
	my $genome = shift;
	return unless $genome;

	my $owner_list = $genome->owner_list;
	return unless ($owner_list);

	my $owner_group = $owner_list->group;
	return unless ($owner_group);

	my $creator = $owner_group->creator;
	return unless ($creator);

	return $creator;
}

sub get_sequence_types {
	my $type_id = shift;
	
	my $html;
	foreach my $type ( $coge->resultset('GenomicSequenceType')->all() ) {
		$html .= '<option value="' . $type->id . '"' . (defined $type_id && $type_id == $type->id ? ' selected' : '') . '>' . $type->info . '</option>';
	}
	
	return $html;
}

sub get_sources {
	#my %opts = @_;
	
	my %unique;
	foreach ($coge->resultset('DataSource')->all()) {
		$unique{$_->name}++;
	}
	
	return encode_json([sort keys %unique]);
}

sub get_groups_with_access {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my @groups = map { $_->group } $genome->lists;
	
	my @rows;
	foreach my $group (sort {$a->name cmp $b->name} @groups) {
		#next if (!$USER->is_admin && !$is_user && !$group->has_member($USER));
		next if ($group->is_owner and !$USER->is_admin);

		my $id = $group->id;
		my %row;
		$row{GROUP_NAME} = qq{<span class="link" onclick='window.open("GroupView.pl?ugid=$id")'>} . $group->name . "</span>";
		$row{GROUP_ROLE} = '(' . $group->role->name . ')';
		$row{GROUP_DESC} = $group->description if $group->description;
		#$row{GROUP_ID} = $id;
		push @rows, \%row;
	}
	push @rows, {GROUP_NAME => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_GROUPS => 1, GROUP_LOOP => \@rows );
	return $template->output;
}

sub get_users_with_access {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless $genome;
	}

	my @rows;
	my $owner = get_genome_owner($genome);
	if ($owner) {
		my %users;
		foreach my $list ($genome->lists) {
			foreach my $user ($list->group->users) {
				$users{$user->id} = $user;
			}
		}
		
		
		foreach my $u (sort {$a->display_name cmp $b->display_name} values %users) {
			my $id = $u->id;
			my $is_owner = ($id == $owner->id);

			my %row;
			$row{USER_INFO} = qq{<span class="link" onclick='window.open("User.pl?uid=$id")'>} . $u->display_name . "</span>" . ($is_owner ? ' (owner)' : '');
			$row{USER_ID} = $id if (not $is_owner);
			push @rows, \%row;
		}
	}
	push @rows, {USER_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_USERS => 1, USER_LOOP => \@rows );
	return $template->output;
}

sub get_available_users {
	my %opts = @_;
	my $gid = $opts{gid};
	return unless $gid;
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my %exists;
	foreach my $list ($genome->lists) {
		foreach my $user ($list->group->users) {
			$exists{$user->id} = $user;
		}
	}

	my @all_users;
	foreach ($coge->resultset('User')->all) {
		next if $_->id == $USER->id; # skip self
		push @all_users, { UID => $_->id, UNAME => $_->display_name } if (!$exists{$_->id});
	}	
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( ADD_USERS => 1 );
	$template->param( GID         => $genome->id );
	$template->param( ALL_UID_LOOP => [sort {$a->{UNAME} cmp $b->{UNAME}} @all_users] );

	return $template->output;
}

sub get_available_groups {
	my %opts = @_;
	my $gid = $opts{gid};
	return unless $gid;
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my %exists = map { $_->group->id => 1 } $genome->lists;

	my @all_groups;
	foreach ($USER->groups) {
		#next if $_->id == $USER->id; # skip self
		push @all_groups, { UGID => $_->id, UGNAME => $_->info } if (!$exists{$_->id});
	}
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( ADD_GROUPS => 1 );
	$template->param( GID        => $genome->id );
	$template->param( ALL_UGID_LOOP => [sort {$a->{UGNAME} cmp $b->{UGNAME}} @all_groups] );

	return $template->output;
}

sub create_shared_group_and_list {
	my %opts = @_;
	my $dbx = $opts{dbx};
	my $uid = $opts{uid};

	my $user = $coge->resultset('User')->find($uid);
	return unless ($user);

	my $shared_group = $coge->resultset('UserGroup')->find_or_create({
		creator_user_id => 0,
		name => $user->name,
		description => "Shared group",
		role_id => 4, #FIXME hardcoded value for "reader" role
		locked => 1
	});
	return unless $shared_group;

	my $conn = $coge->resultset('UserGroupConnector')->find_or_create({
		user_id => $uid,
		user_group_id => $shared_group->id
	});
	return unless $conn;

	my $shared_list = $coge->resultset('List')->find_or_create({
		name => $user->name,
		description => "Shared list",
		list_type_id => 7, #FIXME hardcoded value for "shared" list type
		user_group_id => $shared_group->id,
		locked => 1,
		restricted => 1
	});
	return unless $shared_list;

	return $shared_group;
}

sub add_user {
	my %opts = @_;
	my $gid = $opts{gid};
	my $uid = $opts{uid};
	print STDERR "add_user: $gid $uid\n";
	return unless ($gid and $uid);
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless $genome;

	my $user = $coge->resultset('User')->find($uid);
	return unless $user;

	my $shared_group = $USER->shared_group;
	if (!$shared_group) {
		$shared_group = create_shared_group_and_list(dbx => $coge, uid => $USER->id);
		return unless $shared_group;
	}

	my $conn = $coge->resultset('ListConnector')->find_or_create({
		parent_id => $shared_group->shared_list->id,
		child_id => $genome->id,
		child_type => 2 #FIXME hardcoded "genome" child type
	});
	return unless $conn;

	$conn = $coge->resultset('UserGroupConnector')->find_or_create({
		user_id => $uid,
		user_group_id => $shared_group->id
	});
	return unless $conn;

	return 1;
}

sub remove_user {
	my %opts = @_;
	my $gid = $opts{gid};
	my $uid = $opts{uid};
	return unless ($gid and $uid);
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my $user = $coge->resultset('User')->find($uid);
	return unless $user;

	my $shared_group = $USER->shared_group;
	return unless $shared_group;

	my $conn = $coge->resultset('UserGroupConnector')->find({
		user_id => $uid,
		user_group_id => $shared_group->id
	});
	return unless $conn;
	$conn->delete;

	return 1;
}

sub add_to_notebook {
	my %opts = @_;
	my $gid = $opts{gid};
	return unless $gid;
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( ADD_TO_NOTEBOOK => 1 );
	$template->param( GID      => $genome->id );

	return $template->output;
}

sub add_genome_to_notebook {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $nid  = $opts{nid};
	# print STDERR "add_genome_to_notebook: $gid $nid\n";
	return unless ($gid and $nid);
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my $list = $coge->resultset('List')->find($nid);
	return unless ($list);
	
	my $conn = $coge->resultset('ListConnector')->find_or_create({
		parent_id => $nid,
		child_id => $gid,
		child_type => 2, #FIXME hardcoded "genome" child type
	});
	return unless $conn;
	
	return 1;
}

sub remove_genome_from_notebook {
	my %opts  = @_;
	my $gid  = $opts{gid};
	my $nid  = $opts{nid};
	# print STDERR "remove_genome_from_notebook: $gid $nid\n";
	return unless ($gid and $nid);

	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	my $list = $coge->resultset('List')->find($nid);
	return unless ($list);
	
	my $conn = $coge->resultset('ListConnector')->find_or_create({
		parent_id => $nid,
		child_id => $gid,
		child_type => 2, #FIXME hardcoded "genome" child type
	});
	return unless $conn;
	$conn->delete;
	
	return 1;
}

sub search_notebooks { 
	my %opts = @_;
	my $gid 		= $opts{gid};
	my $search_term	= $opts{search_term};
	my $timestamp 	= $opts{timestamp};
#	print STDERR "$ugid $search_term $timestamp\n";
	return unless $gid;
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless ($genome);

	# Get lists already containing this genome
	my %exists = map { $_->id => 1 } $genome->lists;
	
	my @lists;
	my $num_results;
	my $group_str = join(',', map { $_->id } $USER->groups);
	
	# Try to get all items if blank search term
	if (not $search_term) {
		# Get all lists
		if ($USER->is_admin) {
			$num_results = $coge->resultset("List")->count;
			if ($num_results < $MAX_SEARCH_RESULTS) {
				@lists = $coge->resultset("List")->all;
			}
		}
		else {
			my $sql = "locked=0 AND (restricted=0 OR user_group_id IN ( $group_str ))";
			$num_results = $coge->resultset("List")->count_literal($sql);
			if ($num_results < $MAX_SEARCH_RESULTS) {
				@lists = $coge->resultset("List")->search_literal($sql);
			}
		}
	}
	# Perform search
	else {
		$search_term = '%'.$search_term.'%';
		if ($USER->is_admin) {
			@lists = $coge->resultset("List")->search(
				\[ 'name LIKE ? OR description LIKE ?', 
				['name', $search_term ], ['description', $search_term] ]);		
		}
		else {
			@lists = $coge->resultset("List")->search_literal("locked=0 AND (restricted=0 OR user_group_id IN ( $group_str )) AND (name LIKE '$search_term' OR description LIKE '$search_term')");
		}
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
	my $html = '';
	foreach my $l (sort listcmp @lists) {
		my $disable = $exists{$l->id} ? "disabled='disabled'" : '';
		$html .= "<option $disable value='" . $l->id . "'>" . $l->info . "</option><br>\n";	
	}
	
	return encode_json({timestamp => $timestamp, html => $html});
}

sub get_notebook_preview {
	my %opts  = @_;
	my $nid = $opts{nid};

	my $list = $coge->resultset('List')->find($nid);
	my $html = '';
	if ($list) {
		my $contents_summary = $list->contents_summary_html;
		$html .= "Notebook <b>'" . $list->name . "'</b> (id" . $list->id . ') ';
		if ($contents_summary) {
			$html .= 'contains:<div style="padding-left: 10px; padding-top: 3px;">' . $contents_summary . '</div>';
		}
		else {
			$html .= 'is empty.';
		}	
	}

	return $html;	
}

sub get_notebooks {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my @lists = $genome->lists;

	my @rows;
	foreach my $list (sort listcmp @lists) {
		#next if ($list->restricted && !$USER->is_admin && !$is_user && !$USER->has_access(list => $list));
		next if ($list->is_owner and !$USER->is_admin);

		my $id = $list->id;
		my %row;
		$row{NOTEBOOK_INFO} = qq{<span class="link" onclick='window.open("NotebookView.pl?lid=$id")'>} . $list->info . "</span>";
		$row{NOTEBOOK_ID} = $id;

		push @rows, \%row;
	}
	push @rows, {NOTEBOOK_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_NOTEBOOKS => 1, NOTEBOOK_LOOP => \@rows );
	return $template->output;
}

sub get_experiments {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my @experiments = $genome->experiments;

	my @rows;
	foreach my $exp (sort experimentcmp @experiments) {
		next if ($exp->deleted);
		#next if ($exp->restricted && !$USER->is_admin && !$is_user && !$USER->has_access(experiment => $exp));

		my $id = $exp->id;
		my %row;
		$row{EXPERIMENT_INFO} = qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=$id")'>} . $exp->info . "</span>";
		
		push @rows, \%row;
	}
	push @rows, {EXPERIMENT_INFO => '<span style="font-style:italic;color:gray;">None</span>'} unless (@rows);

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_EXPERIMENTS => 1, EXPERIMENT_LOOP => \@rows );
	return $template->output;
}

sub generate_html {
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= ' ' . $USER->last_name if ( $USER->first_name && $USER->last_name );

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( 
		PAGE_TITLE 	=> $PAGE_TITLE,
		HELP       	=> '/wiki/index.php?title=' . $PAGE_TITLE . '.pl',
		USER     	=> $name,
		LOGO_PNG 	=> $PAGE_TITLE . "-logo.png",
		DATE     	=> $DATE,
		BODY 		=> generate_body(),
		ADJUST_BOX 	=> 1 
	);	
	$template->param( LOGON => 1 ) unless ($USER->user_name eq "public");

	return $template->output;
}

sub generate_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( MAIN => 1, PAGE_NAME => $PAGE_TITLE . '.pl' );
	
	my $gid = $FORM->param('gid');
	return "No genome specified" unless $gid;

	my $genome = $coge->resultset('Genome')->find($gid);
	return "Genome id$gid not found" unless ($genome);	

	return "Access denied" unless (!$genome->restricted or $USER->is_admin or $USER->has_access_to_genome($genome));

	my ($first_chr) = sort {$b->sequence_length <=> $a->sequence_length} $genome->chromosomes;

	$template->param(
		GID 			=> $gid,
		CHR 			=> ($first_chr ? $first_chr : ''),
		GENOME_INFO 	=> get_genome_info(genome => $genome),
	 	GENOME_DATA 	=> get_genome_data(genome => $genome),
		#GROUPS 			=> get_groups_with_access(genome => $genome),
		#USERS 			=> get_users_with_access(genome => $genome),
		#NOTEBOOKS 		=> get_notebooks(genome => $genome),
		EXPERIMENTS 	=> get_experiments(genome => $genome) 
	);

	$template->param( ADMIN_AREA => 1 ) if ($USER->is_admin);

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

# FIXME this routine is duplicated elsewhere
sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}
