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

my $node_types = CoGeX::node_types();

$P         = CoGe::Accessory::Web::get_defaults("$ENV{HOME}/coge.conf");
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
	edit_genome_info		=> \&edit_genome_info,
	update_genome_info		=> \&update_genome_info,
	update_owner			=> \&update_owner,
	search_organisms		=> \&search_organisms,
	search_users			=> \&search_users,
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
	my $genome = $opts{genome};
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
		TYPE => $genome->type->info,
		SOURCE => join(',', map { $_->name } $genome->source),
		RESTRICTED => ($genome->restricted ? 'Yes' : 'No'),
		NAME => $genome->name,
		DESCRIPTION => $genome->description,
		DELETED => $genome->deleted
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

sub update_owner {
	my %opts = @_;
	my $gid = $opts{gid};
	my $user_name = $opts{user_name};
	return unless $gid and $user_name;
	
	# Admin-only function
	return unless $USER->is_admin;

	# Make new user owner of genome
	my $user = $coge->resultset('User')->find( { user_name => $user_name } );
	unless ($user) {
		return "error finding user '$user_name'\n";
	}
	my $conn = $coge->resultset('UserConnector')->find_or_create(
	  { parent_id => $user->id,
		parent_type => $node_types->{user},
		child_id => $gid,
		child_type => $node_types->{genome},
		role_id => 2 # FIXME hardcoded
	  } );
	unless ($conn) {
		return "error creating user connector\n";
	}
	
	# Remove admin user as owner
	$conn = $coge->resultset('UserConnector')->find(
	  { parent_id => $USER->id,
		parent_type => $node_types->{user},
		child_id => $gid,
		child_type => $node_types->{genome},
		role_id => 2 # FIXME hardcoded
	  } );
	if ($conn) {
		$conn->delete;
	}
	
	return;
}

sub get_genome_data {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);
	$gid = $genome->id if $genome;
	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	#ERIC added these 1/31/2013
	#TODO: this should happen in the genome object
	my $seq_file = $genome->file_path;
	my $cogedir = $P->{COGEDIR};
	my $cogeurl = $P->{URL};
	$seq_file =~ s/$cogedir/$cogeurl/i;
	#TODO: this should happen in the genome object
	my $download .= qq{<a class=link href='$seq_file' target="_new">Fasta Sequences</a>};
	$download .= qq{&nbsp|&nbsp};
	$download .= qq{<span class=link onclick="\$('#gff_export').dialog('option', 'width', 400).dialog('open')">Export GFF</span>};
	$download .= qq{&nbsp|&nbsp};
	$download .= qq{<span class=link onclick="export_tbl('$gid')"">Export TBL</span>};
	$download .= qq{&nbsp|&nbsp};
	$download .= qq{<span class=link onclick="export_bed('$gid')"">Export bed</span>};	
	#TODO: this should happen in the genome object
	my $links = "<a href='OrganismView.pl?dsgid=$gid' target=_new>OrganismView</a>";
	$links .= qq{&nbsp|&nbsp};
	$links .= "<a href='CodeOn.pl?dsgid=$gid' target=_new>CodeOn</a>";
	$links .= qq{&nbsp|&nbsp};
	$links .= qq{<span class='link' onclick="window.open('GenomeInfo.pl?gid=$gid');">GenomeInfo</span>};
	$links .= qq{&nbsp|&nbsp};
	$links .= qq{<span class='link' onclick="window.open('SynMap.pl?dsgid1=$gid;dsgid2=$gid');">SynMap</span>};
	$links .= qq{&nbsp|&nbsp};
	$links .= qq{<span class='link' onclick="window.open('CoGeBlast.pl?dsgid=$gid');">CoGeBlast</span>};
	$template->param(
		DO_GENOME_DATA => 1,
		CHROMOSOME_COUNT => commify($genome->chromosome_count()),
		LENGTH => commify($genome->length),
		GID => $genome->id,
		DOWNLOAD=>$download,
		LINKS=>$links,
	);


	return $template->output;
}

sub get_genome_download_links
  {
    my $genome = shift;
    
  }

sub get_sequence_types {
	my $type_id = shift;
	
	my $html;
	foreach my $type ( sort {$a->info cmp $b->info} $coge->resultset('GenomicSequenceType')->all() ) {
		$html .= '<option value="' . $type->id . '"' 
			  . (defined $type_id && $type_id == $type->id ? ' selected' : '') . '>' 
			  . $type->info 
			  . '</option>';
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
	return "" unless @experiments;

	my @rows;
	foreach my $exp (sort experimentcmp @experiments) {
		next if ($exp->deleted);
		#next if ($exp->restricted && !$USER->is_admin && !$is_user && !$USER->has_access(experiment => $exp));

		my $id = $exp->id;
		my %row;
		$row{EXPERIMENT_INFO} = qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=$id")'>} . $exp->info . "</span>";
		
		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_EXPERIMENTS => 1,
					  EXPERIMENT_LOOP => \@rows
	);
	return $template->output;
}

sub get_datasets {
	my %opts = @_;
	my $gid  = $opts{gid};
	my $genome  = $opts{genome};
	return unless ($gid or $genome);

	unless ($genome) {
		$genome = $coge->resultset('Genome')->find($gid);
		return unless ($genome);
	}

	my @datasets = $genome->datasets;
	return "" unless @datasets;

	my @rows;
	foreach my $ds (sort @datasets) {
		my %row;
		$row{DATASET_INFO} = qq{<span>} . $ds->info . "</span>";
		push @rows, \%row;
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
	$template->param( DO_DATASETS => 1,
					  DATASET_LOOP => \@rows
	);
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

	my ($first_chr) = $genome->chromosomes;

	$template->param(
		GID 			=> $gid,
		CHR 			=> ($first_chr ? $first_chr : ''),
		GENOME_INFO 	=> get_genome_info(genome => $genome),
	 	GENOME_DATA 	=> get_genome_data(genome => $genome),
		EXPERIMENTS 	=> get_experiments(genome => $genome) 
	);

	if ($USER->is_admin) {
		$template->param( ADMIN_AREA => 1,
						  DATASETS => get_datasets(genome => $genome) );
	}

	return $template->output;
}

# FIXME this routine is duplicated elsewhere
sub experimentcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->name cmp $b->name || $b->id cmp $a->id
}

# FIXME this routine is duplicated elsewhere
sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}
