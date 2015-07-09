#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadGenome.pl, replicate changes when applicable.

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Accessory::IRODS;
use CoGe::Core::Storage qw(create_annotation_dataset get_workflow_paths get_irods_path get_irods_file  get_upload_path);
use CoGe::Core::Genome qw(genomecmp);
use HTML::Template;
use JSON::XS;
use URI::Escape::JavaScript qw(escape unescape);
use File::Path;
use File::Copy;
use File::Basename;
use File::Slurp;
use File::Spec::Functions qw(catdir catfile);
use File::Listing qw(parse_dir);
use LWP::Simple;
use URI;
use Sort::Versions;
no warnings 'redefine';

use vars qw(
  $P $PAGE_TITLE $LINK $EMBED $TEMPDIR $USER $coge $FORM
  %FUNCTION $MAX_SEARCH_RESULTS $CONFIGFILE $LOAD_ID $WORKFLOW_ID
);

$PAGE_TITLE = 'LoadAnnotation';

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$CONFIGFILE = $ENV{COGE_HOME} . '/coge.conf';

# Get workflow_id and load_id for previous load if specified.  Otherwise
# generate a new load_id for data upload.
$WORKFLOW_ID = $FORM->Vars->{'wid'} || $FORM->Vars->{'job_id'}; # wid is new name, job_id is legacy name
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = get_upload_path($USER->name, $LOAD_ID);

$EMBED = $FORM->param('embed');

$MAX_SEARCH_RESULTS = 100;

%FUNCTION = (
    irods_get_path  => \&irods_get_path,
    irods_get_file  => \&irods_get_file,
    load_from_ftp   => \&load_from_ftp,
    ftp_get_file    => \&ftp_get_file,
    upload_file     => \&upload_file,
    load_annotation => \&load_annotation,
    get_sources     => \&get_sources,
    create_source   => \&create_source,
    search_genomes  => \&search_genomes,
    get_load_log    => \&get_load_log,
	check_login     => \&check_login,
	send_error_report => \&send_error_report
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub generate_html {
    $EMBED = $FORM->param('embed');
    
    # Check for finished result # mdb removed 3/4/15 no longer auto-redirect, make user select result
#    if ($JOB_ID) {
#        my $log = get_load_log(workflow_id => $JOB_ID);
#        if ($log) {
#            my $res = decode_json($log);
#            if ($res->{genome_id}) {
#                my $url = 'GenomeInfo.pl?embed=' . $EMBED . '&gid=' . $res->{genome_id};
#                print $FORM->redirect(-url => $url);
#            }
#        }
#    }
    
    my $template;

    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param( PAGE_TITLE => $PAGE_TITLE,
					      TITLE      => "LoadAnnotation",
        				  PAGE_LINK  => $LINK,
        				  HOME       => $P->{SERVER},
                          HELP       => 'LoadAnnotation',
                          WIKI_URL   => $P->{WIKI_URL} || '',
                          USER       => $USER->display_name || '' );
        $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
        
        my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
        $link = CoGe::Accessory::Web::get_tiny_link( url => $link );
    
        $template->param( ADJUST_BOX => 1 );
        $template->param( ADMIN_ONLY => $USER->is_admin );
        $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    }
    
    $template->param( BODY => generate_body() );
    return $template->output;
}

sub generate_body {
    if ( $USER->user_name eq 'public' ) {
        my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
        $template->param( LOGIN     => 1 );
        return $template->output;
    }

    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );

    my $gid;
    $gid = $FORM->param('gid') if defined $FORM->param('gid');
    if ($gid) {
        my $genome = $coge->resultset('Genome')->find($gid);
        #TODO check permissions
        if ($genome) {
            $template->param(
                GENOME_NAME => $genome->info,
                GENOME_ID   => $genome->id
            );
        }
    }

    $template->param(
        MAIN          => 1,
        PAGE_TITLE    => $PAGE_TITLE,
        EMBED         => $EMBED,
        LOAD_ID       => $LOAD_ID,
        WORKFLOW_ID   => $WORKFLOW_ID,
        API_BASE_URL  => 'api/v1/', #TODO move into config file or module
        HELP_URL      => 'https://genomevolution.org/wiki/index.php/LoadAnnotation',
        SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
        DEFAULT_TAB              => 0,
        FILE_SELECT_SINGLE       => 1,
        DISABLE_IRODS_GET_ALL    => 1,
        MAX_IRODS_LIST_FILES     => 1000,
        MAX_IRODS_TRANSFER_FILES => 30,
        MAX_FTP_FILES            => 30,
        USER                     => $USER->user_name,
    );
    $template->param( SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
                      SPLASH_CONTENTS    => 'This page allows you to load genome sequences in FASTA file format.' );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

sub irods_get_path {
    my %opts = @_;
    my $path = $opts{path};
    $path = unescape($path);
    #print STDERR "irods_get_path ", $path, "\n";
    
    unless ($path) {
        my $username = $USER->name;
        my $basepath = CoGe::Accessory::Web::get_defaults()->{IRODSDIR};
        $basepath =~ s/\<USER\>/$username/;
        $path = $basepath;
    }
    
    my $result = get_irods_path($path);
    
    if ($result->{error}) {
        # Test for recent new account.  The iPlant IRODS isn't ready for a few
        # minues the first time a user logs into CoGe.
        # mdb added 3/31/14
        my $isNewAccount = 0;
        if ($USER->date ne '0000-00-00 00:00:00') {
            my $dt_user = DateTime::Format::MySQL->parse_datetime( $USER->date );
            my $dt_now = DateTime->now( time_zone => 'America/Phoenix' );
            my $diff = $dt_now->subtract_datetime($dt_user);
            my ( $years, $months, $days, $hours, $minutes ) = $diff->in_units('years', 'months', 'days', 'hours', 'minutes');
            $isNewAccount = (!$years && !$months && !$days && !$hours && $minutes < 5) ? 1 : 0;
        }

        # Send support email
        if (!$isNewAccount) {
            my $email = $P->{SUPPORT_EMAIL};
            my $body =
                "irods ils command failed\n\n"
              . 'User: '
              . $USER->name . ' id='
              . $USER->id . ' '
              . $USER->date . "\n\n"
              . $result->{error} . "\n\n"
              . $P->{SERVER};
            CoGe::Accessory::Web::send_email(
                from    => $email,
                to      => $email,
                subject => "System error notification",
                body    => $body
            );
        }
        return encode_json($result);
    }

    return encode_json({ path => $path, items => $result->{items} });
}

sub irods_get_file {
    my %opts = @_;
    my $path = $opts{path};
    
    $path = unescape($path);
    
    my $result = get_irods_file($path, $TEMPDIR);

    return encode_json( { path => $result->{path}, size => $result->{size} } );
}

sub load_from_ftp {
    my %opts = @_;
    my $url  = $opts{url};

    my @files;

    my ($content_type) = head($url);
    if ($content_type) {
        if ( $content_type eq 'text/ftp-dir-listing' ) {    # directory
            my $listing = get($url);
            my $dir     = parse_dir($listing);
            foreach (@$dir) {
                my ( $filename, $filetype, $filesize, $filetime, $filemode ) =
                  @$_;
                if ( $filetype eq 'f' ) {
                    push @files, { name => $filename, url => $url . $filename };
                }
            }
        }
        else {                                              # file
            my ($filename) = $url =~ /([^\/]+?)(?:\?|$)/;
            push @files, { name => $filename, url => $url };
        }
    }
    else {    # error (url not found)
        return;
    }

    return encode_json( \@files );
}

sub ftp_get_file {
    my %opts      = @_;
    my $url       = $opts{url};
    my $username  = $opts{username};
    my $password  = $opts{password};

    #my ( $type, $filepath, $filename ) = $url =~ /^(ftp|http):\/\/(.+)\/(\S+)$/; # mdb removed 1/6/14, issue 274
	# mdb added 1/6/14, issue 274
	my $uri = URI->new($url);
	my $type = $uri->scheme;
	my ($filename, $filepath) = fileparse($uri->path);
	$filepath = $uri->host . $filepath;

    # print STDERR "$type $filepath $filename $username $password\n";
    return unless ( $type and $filepath and $filename );

    my $path         = 'ftp/' . $filepath . '/' . $filename;
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
    my $request = HTTP::Request->new( GET => $url );
    $request->authorization_basic( $username, $password )
      if ( $username and $password );

    #print STDERR "request uri: " . $request->uri . "\n";
    $request->content_type("text/xml; charset=utf-8");
    my $response = $ua->request($request);
    if ( $response->is_success() ) {

        #my $header = $response->header;
        my $result = $response->content;

        #print STDERR "content: <begin>$result<end>\n";
        open( my $fh, ">$fullfilepath/$filename" );
        if ($fh) {
            binmode $fh;    # could be binary data
            print $fh $result;
            close($fh);
        }
    }
    else {                  # error
        my $status = $response->status_line();
        print STDERR "status_line: $status\n";
        return encode_json(
            {
                path      => $path,
                size      => "Failed: $status"
            }
        );
    }

    return encode_json(
        {
            path      => $path,
            size      => -s $fullfilepath . '/' . $filename
        }
    );
}

sub upload_file {
    my %opts      = @_;
    my $filename  = '' . $FORM->param('input_upload_file');
    my $fh        = $FORM->upload('input_upload_file');

    #   print STDERR "upload_file: $filename\n";

    my $size = 0;
    my $path;
    if ($fh) {
        my $tmpfilename = $FORM->tmpFileName( $FORM->param('input_upload_file') );
        $path = catfile('upload', $filename);
        my $targetpath = catdir($TEMPDIR, 'upload');
        mkpath($targetpath);
        $targetpath = catfile($targetpath, $filename);

        #print STDERR "temp files: $tmpfilename $targetpath\n";
        copy( $tmpfilename, $targetpath );
        $size = -s $fh;
    }

    return encode_json(
        {
            filename  => $filename,
            path      => $path,
            size      => $size
        }
    );
}

sub check_login {
	print STDERR $USER->user_name . ' ' . int($USER->is_public) . "\n";
	return ($USER && !$USER->is_public);
}

sub load_annotation {
    my %opts        = @_;
    my $name        = $opts{name};
    my $description = $opts{description};
    my $link        = $opts{link};
    my $version     = $opts{version};
    my $source_name = $opts{source_name};
    my $user_name   = $opts{user_name};
    my $gid         = $opts{gid};
    my $items       = $opts{items};

    # print STDERR "load_annotation: name=$name description=$description version=$version gid=$gid\n";

    # Check login
    if ( !$user_name || !$USER->is_admin ) {
        $user_name = $USER->user_name;
    }
    if ($user_name eq 'public') {
        return encode_json({ error => 'Not logged in' });
    }
    
    # Check genome
    my $genome = $coge->resultset('Genome')->find($gid);
    return unless $genome;

    # Check data items
    return encode_json({ error => "No data items" }) unless $items;
    $items = decode_json($items);
    #print STDERR Dumper $items;

    # mdb added issue 309 - workaround here because checking perms in search_genomes() is too slow
    return encode_json({ error => "You do not have permission to modify this genome" })
        unless ($USER->is_admin || $USER->is_owner_editor( dsg => $gid ));

    # Setup paths to files
    my @files = map { $TEMPDIR . $_->{path} } @$items;

    # Submit workflow to add genome
    my ($workflow_id, $error_msg) = create_annotation_dataset(
        user => $USER,
        metadata => {
            name => $name,
            description => $description,
            link => $link,
            version => $version,
            source_name => $source_name,
            genome_id => $gid
        },
        files => \@files
    );
    unless ($workflow_id) {
        return encode_json({ error => "Workflow submission failed: " . $error_msg });
    }

	# Get tiny link
    my $tiny_link = CoGe::Accessory::Web::get_tiny_link(
        url => $P->{SERVER} . "$PAGE_TITLE.pl?job_id=" . $workflow_id . "&embed=" . $EMBED
    );
    
    # Log it
    my $info;
    $info .= " v" . $version;
    $info .= ' for genome <i>"' . $genome->organism->name;
    $info .= " (" . $genome->name . ")" if $genome->name;
    $info .= ": " . $genome->description if $genome->description;
    $info .= " (v" . $genome->version . ")";
    $info .= '"</i>';
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        workflow_id => $workflow_id,
        user_id     => $USER->id,
        page        => "LoadAnnotation",
        description => 'Load annotation ' . $info,
        link        => $tiny_link
    );

    return encode_json({ job_id => $workflow_id, link => $tiny_link });
}

sub get_load_log {
    my %opts         = @_;
    my $workflow_id = $opts{workflow_id};
    return unless $workflow_id;
    #TODO authenticate user access to workflow

    my (undef, $results_path) = get_workflow_paths($USER->name, $workflow_id);
    return unless (-r $results_path);

    my $result_file = catfile($results_path, '1');
    return unless (-r $result_file);

    my $result = CoGe::Accessory::TDS::read($result_file);
    return unless $result;

    return encode_json(
        {
            genome_id   => $result->{genome_id},
            dataset_id  => $result->{dataset_id},
        }
    );
}

sub search_genomes {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    my $id = $search_term;
    $search_term = '%' . $search_term . '%';

    # Get all matching organisms
    my @organisms = $coge->resultset("Organism")->search(
        \[
            'name LIKE ? OR description LIKE ?',
            [ 'name',        $search_term ],
            [ 'description', $search_term ]
        ]
    );

    # Get all matching genomes
    my @genomes = $coge->resultset("Genome")->search(
        \[
            'genome_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'genome_id',   $id ],
            [ 'name',        $search_term ],
            [ 'description', $search_term ]
        ]
    );

	# Combine matching genomes with matching organism genomes, preventing duplicates
    my %unique;
    map {
        $unique{ $_->id } = $_ if ($USER->has_access_to_genome($_))
    } @genomes;
    foreach my $organism (@organisms) {
        map {
            $unique{ $_->id } = $_ if ($USER->has_access_to_genome($_))
        } $organism->genomes;
    }

    # Limit number of results displayed
    if ( keys %unique > $MAX_SEARCH_RESULTS ) {
        return encode_json( { timestamp => $timestamp, items => undef } );
    }

    my @items;
    foreach ( sort genomecmp values %unique ) {    #(keys %unique) {
        push @items, { label => $_->info, value => $_->id };
    }

    return encode_json( { timestamp => $timestamp, items => \@items } );
}

sub get_sources {

    #my %opts = @_;

    my %unique;
    foreach ( $coge->resultset('DataSource')->all() ) {
        $unique{ $_->name }++;
    }

    return encode_json( [ sort keys %unique ] );
}

sub create_source {
    my %opts = @_;
    my $name = $opts{name};
    return unless $name;
    my $desc = $opts{desc};
    my $link = $opts{link};
    $link =~ s/^\s+//;
    $link = 'http://' . $link if ( not $link =~ /^(\w+)\:\/\// );

    my $source =
      $coge->resultset('DataSource')
      ->find_or_create(
        { name => $name, description => $desc, link => $link } );
    return unless ($source);

    return $name;
}

sub get_debug_log {
    my %opts         = @_;
    my $workflow_id = $opts{workflow_id};
    return unless $workflow_id;
    #TODO authenticate user access to workflow

    my (undef, $results_path) = get_workflow_paths($USER->name, $workflow_id);
    return unless (-r $results_path);

    my $result_file = catfile($results_path, 'debug.log');
    return unless (-r $result_file);

    my $result = read_file($result_file);
    return $result;
}

sub send_error_report {
    my %opts = @_;
    my $load_id = $opts{load_id};
    my $job_id = $opts{job_id};

    # Get the staging directory
    my ($staging_dir, $result_dir) = get_workflow_paths($USER->name, $job_id);

    my $url = $P->{SERVER} . "$PAGE_TITLE.pl?";
    $url .= "job_id=$job_id;" if $job_id;
    $url .= "load_id=$load_id";

    my $email = $P->{SUPPORT_EMAIL};

    my $body =
        "Load failed\n\n"
        . 'For user: '
        . $USER->name . ' id='
        . $USER->id . ' '
        . $USER->date . "\n\n"
        . "staging_directory: $staging_dir\n\n"
        . "result_directory: $result_dir\n\n"
        . "tiny link: $url\n\n";

    $body .= get_debug_log(workflow_id => $job_id);

    CoGe::Accessory::Web::send_email(
        from    => $email,
        to      => $email,
        subject => "Load error notification from $PAGE_TITLE",
        body    => $body
    );
}
