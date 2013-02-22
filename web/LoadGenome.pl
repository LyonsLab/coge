#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadExperiment.pl, replicate changes when applicable.

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
	%FUNCTION $MAX_SEARCH_RESULTS $CONFIGFILE
);

$CONFIGFILE = $ENV{HOME} . 'coge.conf';
$P         = CoGe::Accessory::Web::get_defaults($CONFIGFILE);
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_TITLE = 'LoadGenome';

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

$BINDIR = $P->{BINDIR};

$MAX_SEARCH_RESULTS = 100;

#$SIG{'__WARN__'} = sub { };    # silence warnings

%FUNCTION = (
	generate_html			=> \&generate_html,
	irods_get_path			=> \&irods_get_path,
	irods_get_file			=> \&irods_get_file,
	load_from_ftp			=> \&load_from_ftp,
	ftp_get_file			=> \&ftp_get_file,
	upload_file				=> \&upload_file,
	ncbi_search				=> \&ncbi_search,
	load_genome				=> \&load_genome,
	get_sequence_types		=> \&get_sequence_types,
	create_sequence_type	=> \&create_sequence_type,
	create_source			=> \&create_source,
	create_organism			=> \&create_organism,
	search_organisms		=> \&search_organisms,
	search_users			=> \&search_users,
	get_sources				=> \&get_sources,
	get_load_genome_log		=> \&get_load_genome_log,
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

sub irods_get_path {
	my %opts = @_;
	my $path = $opts{path};
	my $timestamp = $opts{timestamp};
	
	my $username = $USER->name;
	my $basepath = $P->{IRODSDIR};
	$basepath =~ s/\<USER\>/$username/;
	$path = $basepath unless $path;

	if ($path !~ /^$basepath/) {
		print STDERR "Attempt to access '$path' denied (basepath='$basepath')\n";
		return;
	}
	
	my $result = CoGe::Accessory::Web::irods_ils($path);
	my $error = $result->{error};
	if ($error) {
		my $email = $P->{SUPPORT_EMAIL};
		CoGe::Accessory::Web::send_email(from => $email, to => $email, subject => "System error notification from $PAGE_TITLE", body => $error);
		return encode_json( { timestamp => $timestamp, error => $error } );
	}
	return encode_json( { timestamp => $timestamp, path => $path, items => $result->{items} } );
}

sub irods_get_file {
	my %opts = @_;
	my $path = $opts{path};
	
	my ($filename) = $path =~ /([^\/]+)\s*$/;
	my ($remotepath) = $path =~ /(.*)$filename$/;
#	print STDERR "get_file $path $filename\n";
	
	my $localpath = 'irods/' . $remotepath;
	my $localfullpath = $TEMPDIR . $localpath;
	$localpath .= '/' . $filename;
	my $localfilepath = $localfullpath . '/' . $filename;
	
	my $do_get = 1;
#	if (-e $localfilepath) {
#		my $remote_chksum = irods_chksum($path);
#		my $local_chksum = md5sum($localfilepath);
#		$do_get = 0 if ($remote_chksum eq $local_chksum);
#		print STDERR "$remote_chksum $local_chksum\n";
#	}
	
	if ($do_get) {
		mkpath($localfullpath);
		CoGe::Accessory::Web::irods_iget($path, $localfullpath);
	}

	return encode_json( { path => $localpath, size => -s $localfilepath } );
}

sub load_from_ftp {
	my %opts = @_;
	my $url = $opts{url};
	
	my @files;
	
	my ($content_type) = head($url);
	if ($content_type) {
		if ($content_type eq 'text/ftp-dir-listing') { # directory
			my $listing = get($url);
			my $dir = parse_dir($listing);
			foreach (@$dir) {
				my ($filename, $filetype, $filesize, $filetime, $filemode) = @$_;
				if ($filetype eq 'f') {
					push @files, { name => $filename, url => $url . $filename };
				}
			}
		}
		else { # file
			my ($filename) = $url =~ /([^\/]+)\s*$/;
			push @files, { name => $filename, url => $url};
		}
	}
	else { # error (url not found)
		return;
	}
		
	return encode_json( \@files );
}

sub ftp_get_file {
	my %opts = @_;
	my $url = $opts{url};
	my $username = $opts{username};
	my $password = $opts{password};
	my $timestamp = $opts{timestamp};
	
	my ($type, $filepath, $filename) = $url =~ /^(ftp|http):\/\/(.+)\/(\S+)$/;
	# print STDERR "$type $filepath $filename $username $password\n";
	return unless ($type and $filepath and $filename);
	
	my $path = 'ftp/' . $filepath . '/' . $filename;
	my $fullfilepath = $TEMPDIR . 'ftp/' . $filepath;
	mkpath($fullfilepath);

	# Simplest method (but doesn't allow login)
#	print STDERR "getstore: $url\n";
#	my $res_code = getstore($url, $fullfilepath . '/' . $filename);
#	print STDERR "response: $res_code\n";
	# TODO check response code here

	# Alternate method with progress callback
#	my $ua = new LWP::UserAgent;
#	my $expected_length;
#	my $bytes_received = 0;
#	$ua->request(HTTP::Request->new('GET', $url),
#		sub {
#			my($chunk, $res) = @_;
#			print STDERR "matt: " . $res->header("Content_Type") . "\n";
#			
#			$bytes_received += length($chunk);
#			unless (defined $expected_length) {
#				$expected_length = $res->content_length || 0;
#			}
#			if ($expected_length) {
#				printf STDERR "%d%% - ",
#				100 * $bytes_received / $expected_length;
#			}
#			print STDERR "$bytes_received bytes received\n";
#				
#			# XXX Should really do something with the chunk itself
##			print STDERR $chunk;
#		});

	# Current method (allows optional login)
	my $ua = new LWP::UserAgent;
	my $request = HTTP::Request->new(GET => $url);
	$request->authorization_basic($username, $password) if ($username and $password);
	#print STDERR "request uri: " . $request->uri . "\n";
	$request->content_type("text/xml; charset=utf-8");
	my $response = $ua->request($request);
	if ($response->is_success()) {
		#my $header = $response->header;
		my $result = $response->content;
		#print STDERR "content: <begin>$result<end>\n";
		open(my $fh, ">$fullfilepath/$filename");
		if ($fh) {
			binmode $fh; # could be binary data
			print $fh $result;
			close($fh);
		}
	}
	else { # error
		my $status = $response->status_line();
		print STDERR "status_line: $status\n";
		return encode_json({ timestamp => $timestamp, path=> $path, size => "Failed: $status" });
	}
	
	return encode_json({ timestamp => $timestamp, path=> $path, size => -s $fullfilepath . '/' . $filename });
}

sub ncbi_search {
	my %opts = @_;
	my $accn = $opts{accn};
	my $timestamp = $opts{timestamp};
	my $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$accn";
	my $result = get( $esearch );
	#print STDERR $result;
	
	my $record = XMLin($result);
	#print STDERR Dumper $record;
	
	my $id = $record->{IdList}->{Id};
	print STDERR "id = $id\n";
	
	my $title;
	if ($id) {
		$esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=$id";
		my $result = get( $esearch );
		#print STDERR $result;
		
		$record = XMLin($result);
		#print STDERR Dumper $record;
		
		foreach (@{$record->{DocSum}->{Item}}) { #FIXME use grep here instead
			if ($_->{Name} eq 'Title') {
				$title = $_->{content};
				print STDERR "title=$title\n";
				last;
			}
		}
	}
	
	return unless $id and $title;
	return encode_json( { timestamp => $timestamp, name => $title, id => $id } );
}

sub upload_file {
	my %opts = @_;
	my $timestamp = $opts{timestamp};
	my $filename = '' . $FORM->param('input_upload_file');
	my $fh = $FORM->upload('input_upload_file');
#	print STDERR "upload_file: $filename\n";
	
	my $size = 0;
	my $path;
	if ($fh) {
		my $tmpfilename = $FORM->tmpFileName( $FORM->param('input_upload_file') );
		$path = 'upload/' . $filename;
		my $targetpath = $TEMPDIR . 'upload/';
		mkpath($targetpath);
		$targetpath .= $filename;
#		print STDERR "temp files: $tmpfilename $targetpath\n";
		copy($tmpfilename, $targetpath);
		$size = -s $fh;
	}

	return encode_json( { timestamp => $timestamp, filename => $filename, path => $path, size => $size } );	
}

sub load_genome {
	my %opts = @_;
	my $name = $opts{name};
	my $description = $opts{description};
	my $link = $opts{link};
	my $version = $opts{version};
	my $type_id = $opts{type_id};
	my $source_name = $opts{source_name};
	my $restricted = $opts{restricted};
	my $organism_id = $opts{organism_id};
	my $user_name = $opts{user_name};
	my $items = $opts{items};
	return unless $items;
	print STDERR "load_genome: organism_id=$organism_id name=$name description=$description version=$version type_id=$type_id restricted=$restricted\n";
	
	$items = decode_json($items);
#	print STDERR Dumper $items;

	if (!$user_name || !$USER->is_admin) {
		$user_name = $USER->user_name;
	}

	# Setup staging area and log file
	my $stagepath = $TEMPDIR . 'staging/';
	my $i;
	for ($i = 1;  -e "$stagepath$i";  $i++) {};
	$stagepath .= $i;
	mkpath $stagepath;
	
	my $logfile = $stagepath . '/log.txt';
	open(my $log, ">$logfile") or die "Error creating log file";
	print $log "Starting load genome $stagepath\n" .
			  "name=$name description=$description version=$version type_id=$type_id restricted=$restricted org_id=$organism_id\n";

	# Verify and decompress files
	my @files;
	foreach my $item (@$items) {
		my $fullpath = $TEMPDIR . $item;
		die "File doesn't exist! $fullpath" if (not -e $fullpath);
		my ($path, $filename) = $item =~ /^(.+)\/([^\/]+)$/;
		my ($fileext) = $filename =~ /\.([^\.]+)$/;
#		print STDERR "$path $filename $fileext\n";
		if ($fileext eq 'gz') { 
			my $cmd = $P->{GUNZIP} . ' ' . $fullpath;
#			print STDERR "$cmd\n";
			`$cmd`;
			$fullpath =~ s/\.$fileext$//;
		}
		#TODO support detecting/untarring tar files also
		push @files, $fullpath;
	}

	print $log "Calling bin/load_genome.pl ...\n";
	my $datadir = $P->{DATADIR} . '/genomic_sequence/';
	my $cmd = "$BINDIR/load_genome.pl " .
			  "-user_name " . $user_name . ' ' .
			  '-name "' . escape($name) . '" ' .
			  '-desc "' . escape($description) . '" ' .
			  '-link "' . escape($link) . '" ' .
			  '-version "' . escape($version) . '" ' .
			  "-type_id $type_id " .
			  "-restricted " . ($restricted eq 'true') . ' ' .
			  "-organism_id $organism_id " .
			  '-source_name "' . escape($source_name) . '" ' .
			  "-staging_dir $stagepath " .
			  "-install_dir $datadir " .
			  '-fasta_files "' . escape(join(',', @files)) . '" ' .
			  "-config $CONFIGFILE";
			  #"-host $DBHOST -port $DBPORT -database $DBNAME -user $DBUSER -password $DBPASS";
	print STDERR "$cmd\n";
	print $log "$cmd\n";	
	close($log);

	if (!defined(my $child_pid = fork())) {
	    die "cannot fork: $!";
	} 
	elsif ($child_pid == 0) {
		print STDERR "child running: $cmd\n";
		`$cmd`;
	    exit;
	} 
	
	return $i;
}

sub get_load_genome_log {
	my %opts = @_;
	my $load_id = $opts{load_id};
	
	my $logfile = $TEMPDIR . "staging/$load_id/log.txt";
	open(my $fh, $logfile) 
		or return encode_json({ status => -1, log => "Error opening log file" });

	my @lines = ();
	my $gid;
	my $status = 0;
	while (<$fh>) {
		push @lines, $1 if ($_ =~ /^log: (.+)/i);
		if ($_ =~ /All done/i) {
			$status = 1;
			last;
		}
		elsif ($_ =~ /log: Added genome id(\d+)/i) {
			$gid = $1;
		}
		elsif ($_ =~ /log: error/i) {
			$status = -1;
			last;
		}
	}
	close($fh);
	
	return encode_json({ status => $status, genome_id => $gid, log => join("<BR>\n", @lines) });
}

sub get_sequence_types {
	my $selected = 1;
	
	my $html;
	foreach my $type ( sort {$a->name cmp $b->name} $coge->resultset('GenomicSequenceType')->all() ) {
		$html .= '<option value="' . $type->id . '"';
		if ($selected && $type->id == $selected) { #$type->name =~ /$selected/i) {
			$html .= ' selected';
			$selected = '';
		}
		$html .= '>' . $type->info . '</option>';
	}
	
	return $html;
}

sub create_sequence_type {
	my %opts = @_;
	my $name = $opts{name};
	my $desc = $opts{desc};
	return unless $name;
	
	my $type = $coge->resultset('GenomicSequenceType')->find_or_create( { name => $name, description => $desc } );
	return unless $type;
	
	return $type->id;
}

sub create_organism {
	my %opts = @_;
	my $name = $opts{name};
	my $desc = $opts{desc};
	return unless $name;
	
	my $organism = $coge->resultset('Organism')->find_or_create( { name => $name, description => $desc } );
	return unless $organism;
	
	return $organism->id;
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
		\[ 'name LIKE ? OR description LIKE ?',  #FIXME security hole: need to check 'restricted'
		['name', $search_term ], ['description', $search_term] ]);

	# Limit number of results displayed
	if (@organisms > $MAX_SEARCH_RESULTS) {
		return encode_json({timestamp => $timestamp, items => undef});
	}
	
	my @results;
	foreach (sort {$a->name cmp $b->name} @organisms) {
		push @results, { 'label' => $_->name, 'value' => $_->id };
	}

	return encode_json({timestamp => $timestamp, items => \@results});
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

sub get_sources {
	#my %opts = @_;
	
	my %unique;
	foreach ($coge->resultset('DataSource')->all()) {
		$unique{$_->name}++;
	}
	
	return encode_json([sort keys %unique]);
}

sub create_source {
	my %opts = @_;
	my $name = $opts{name};
	return unless $name;
	my $desc = $opts{desc};
	my $link = $opts{link};
	$link =~ s/^\s+//;
	$link = 'http://' . $link if (not $link =~ /^(\w+)\:\/\//);	
	
	my $source = $coge->resultset('DataSource')->find_or_create( { name => $name, description => $desc, link => $link } );
	return unless ($source);
	
	return $name;
}

sub generate_html {
	my $html;
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( PAGE_TITLE => $PAGE_TITLE );
	$template->param( HELP       => '/wiki/index.php?title=' . $PAGE_TITLE );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= ' ' . $USER->last_name if ( $USER->first_name && $USER->last_name );
	$template->param( USER     => $name,
					  LOGO_PNG => $PAGE_TITLE . "-logo.png",
					  DATE     => $DATE );
	$template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
	my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
	$link = CoGe::Accessory::Web::get_tiny_link( url => $link );

	$template->param( BODY => generate_body() );
	$template->param( ADJUST_BOX => 1 );

	$html .= $template->output;
	return $html;
}

sub generate_body {
	if ($USER->user_name eq 'public') {
		my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
		$template->param( PAGE_NAME => "$PAGE_TITLE.pl",
						  LOGIN     => 1 );
		return $template->output;
	}
	
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( MAIN => 1,
					  PAGE_NAME => $PAGE_TITLE . '.pl',
					  SUPPORT_EMAIL => $P->{SUPPORT_EMAIL} );

	$template->param( ENABLE_NCBI => 1,
					  DEFAULT_TAB => 0,
					  MAX_IRODS_LIST_FILES => 100,
					  MAX_IRODS_TRANSFER_FILES => 30,
					  MAX_FTP_FILES => 30 );
					  
	$template->param( ADMIN_AREA    => 1 ) if $USER->is_admin;
	
	return $template->output;
}
