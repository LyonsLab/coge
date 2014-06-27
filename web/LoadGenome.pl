#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadExperiment.pl & LoadAnnotation.pl, replicate changes when applicable.

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::IRODS;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage qw(create_genome_from_file create_genome_from_NCBI get_workflow_paths);
use HTML::Template;
use JSON::XS;
use Sort::Versions;
use File::Path qw(mkpath);
use File::Copy qw(copy);
use File::Basename;
use File::Spec::Functions qw( catdir catfile );
use File::Listing qw(parse_dir);
use URI;
use URI::Escape::JavaScript qw(escape);
use LWP::Simple;
use XML::Simple;
use DateTime;
use DateTime::Format::MySQL;
use Data::Dumper;
no warnings 'redefine';

use vars qw(
  $P $PAGE_TITLE $TEMPDIR $user $coge $FORM $LINK $JOB_ID
  %FUNCTION $MAX_SEARCH_RESULTS $CONFIGFILE $LOAD_ID
);

$PAGE_TITLE = 'LoadGenome';

$FORM = new CGI;
( $coge, $user, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$CONFIGFILE = $ENV{COGE_HOME} . '/coge.conf';

$JOB_ID  = $FORM->Vars->{'job_id'};
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = $P->{SECTEMPDIR} . $PAGE_TITLE . '/' . $user->name . '/' . $LOAD_ID . '/';

$MAX_SEARCH_RESULTS = 400;

%FUNCTION = (
    irods_get_path => \&irods_get_path,
    irods_get_file => \&irods_get_file,
    load_from_ftp  => \&load_from_ftp,
    ftp_get_file   => \&ftp_get_file,
    upload_file    => \&upload_file,
    search_ncbi_nucleotide => \&search_ncbi_nucleotide,
    load_genome          => \&load_genome,
    get_sequence_types   => \&get_sequence_types,
    create_sequence_type => \&create_sequence_type,
    create_source        => \&create_source,
    create_organism      => \&create_organism,
    search_organisms     => \&search_organisms,
    search_users         => \&search_users,
    get_sources          => \&get_sources,
    get_load_log         => \&get_load_log,
	check_login			 => \&check_login,
	send_error_report    => \&send_error_report
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub generate_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => $PAGE_TITLE,
    				  PAGE_LINK  => $LINK,
					  HELP       => '/wiki/index.php?title=' . $PAGE_TITLE );
    my $name = $user->user_name;
    $name = $user->first_name if $user->first_name;
    $name .= ' ' . $user->last_name
      if ( $user->first_name && $user->last_name );
    $template->param(
        USER     => $name,
        LOGO_PNG => $PAGE_TITLE . "-logo.png",
    );
    $template->param( LOGON => 1 ) unless $user->user_name eq "public";
    my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
    $link = CoGe::Accessory::Web::get_tiny_link( url => $link );

    $template->param( BODY       => generate_body() );
    $template->param( ADJUST_BOX => 1 );

    $html .= $template->output;
    return $html;
}

sub generate_body {
    if ( $user->user_name eq 'public' ) {
        my $template =
          HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param(
            PAGE_NAME => "$PAGE_TITLE.pl",
            LOGIN     => 1
        );
        return $template->output;
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        MAIN          => 1,
        PAGE_NAME     => $PAGE_TITLE . '.pl',
        LOAD_ID       => $LOAD_ID,
        JOB_ID        => $JOB_ID,
        STATUS_URL    => 'jex/status/',
        SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
        ENABLE_NCBI              => 1,
        DEFAULT_TAB              => 0,
        MAX_IRODS_LIST_FILES     => 100,
        MAX_IRODS_TRANSFER_FILES => 30,
        MAX_FTP_FILES            => 30
    );

    my $oid = $FORM->param("oid");
    my $organism = $coge->resultset('Organism')->find($oid) if $oid;

    if ($organism) {
        $template->param(ORGANISM_NAME => $organism->name);
    }

    $template->param( ADMIN_AREA => 1 ) if $user->is_admin;

    return $template->output;
}

sub irods_get_path {
    my %opts      = @_;
    my $path      = $opts{path};

    my $username = $user->name;
    my $basepath = $P->{IRODSDIR};
    $basepath =~ s/\<USER\>/$username/;
    $path = $basepath unless $path;

    if ( $path !~ /^$basepath/ ) {
        print STDERR "Attempt to access '$path' denied (basepath='$basepath')\n";
        return;
    }

    my $result = CoGe::Accessory::IRODS::irods_ils($path);
    my $error  = $result->{error};

    if ($error) {
        # Test for recent new account.  The iPlant IRODS isn't ready for a few
        # minues the first time a user logs into CoGe.
        # mdb added 3/31/14
        my $isNewAccount = 0;
        if ($user->date ne '0000-00-00 00:00:00') {
            my $dt_user = DateTime::Format::MySQL->parse_datetime( $user->date );
            my $dt_now = DateTime->now( time_zone => 'America/Phoenix' );
            my $diff = $dt_now->subtract_datetime($dt_user);
            my ( $years, $months, $days, $hours, $minutes ) = $diff->in_units('years', 'months', 'days', 'hours', 'minutes');
            $isNewAccount = (!$years && !$months && !$days && !$hours && $minutes < 5) ? 1 : 0;
        }

        if (!$isNewAccount) {
            my $email = $P->{SUPPORT_EMAIL};
            my $body =
                "irods ils command failed\n\n"
              . 'User: '
              . $user->name . ' id='
              . $user->id . ' '
              . $user->date . "\n\n"
              . $error . "\n\n"
              . $P->{SERVER};
            CoGe::Accessory::Web::send_email(
                from    => $email,
                to      => $email,
                subject => "System error notification from $PAGE_TITLE",
                body    => $body
            );
            return encode_json( { error => $error } );
        }
    }
    return encode_json(
        { path => $path, items => $result->{items} } );
}

sub irods_get_file {
    my %opts = @_;
    my $path = $opts{path};

    my ($filename)   = $path =~ /([^\/]+)\s*$/;
    my ($remotepath) = $path =~ /(.*)$filename$/;

    my $localpath     = 'irods/' . $remotepath;
    my $localfullpath = $TEMPDIR . $localpath;
    $localpath .= '/' . $filename;
    my $localfilepath = $localfullpath . '/' . $filename;
    #print STDERR "get_file $path $filename $localfilepath\n";

    my $do_get = 1;

    #	if (-e $localfilepath) {
    #		my $remote_chksum = irods_chksum($path);
    #		my $local_chksum = md5sum($localfilepath);
    #		$do_get = 0 if ($remote_chksum eq $local_chksum);
    #		print STDERR "$remote_chksum $local_chksum\n";
    #	}

    if ($do_get) {
        mkpath($localfullpath);
        CoGe::Accessory::IRODS::irods_iget( $path, $localfullpath );
    }

    return encode_json( { path => $localpath, size => -s $localfilepath } );
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
            my ($filename) = $url =~ /([^\/]+)\s*$/;
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
    my $username  = $opts{username}; # optional
    my $password  = $opts{password}; # optional

    #my ( $type, $filepath, $filename ) = $url =~ /^(ftp|http|https):\/\/(.+)\/(\S+)$/; # mdb removed 1/6/14, issue 274
	# mdb added 1/6/14, issue 274
	my $uri = URI->new($url);
	my $type = $uri->scheme;
	my ($filename, $filepath) = fileparse($uri->path);
	$filepath = $uri->host . $filepath;

    #print STDERR "url=$url type=$type filepath=$filepath filename=$filename ", ($username and $password ? "$username $password" : ''), "\n";
    return unless ( $type and $filepath and $filename );

    my $path         = $type . '/' . $filepath . '/' . $filename;
    my $fullfilepath = $TEMPDIR . '/' . $type . '/' . $filepath;
    #print STDERR "path=$path fullfilepath=$fullfilepath\n";
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
    #           #print STDERR $chunk;
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
            binmode $fh; # could be binary data
            print $fh $result;
            close($fh);
        }
    }
    else { # error
        my $status = $response->status_line();
        #print STDERR "status_line: $status\n";
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

sub search_ncbi_nucleotide { #TODO this can be done client-side instead
    my %opts      = @_;
    my $accn      = $opts{accn};
    my $timestamp = $opts{timestamp};

    my $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$accn";
    my $result = get($esearch);
    #print STDERR $result;
    my $record = XMLin($result);
    #print STDERR Dumper $record;
    my $id = $record->{IdList}->{Id};
    #print STDERR "id = $id\n";

    my $title;
    if ($id) { # GI number
        $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=$id";
        my $result = get($esearch);
        #print STDERR $result;
        $record = XMLin($result);
        #print STDERR Dumper $record;
        foreach ( @{ $record->{DocSum}->{Item} } )
        {    #FIXME use grep here instead
            if ( $_->{Name} eq 'Title' ) {
                $title = $_->{content};
                #print STDERR "title=$title\n";
                last;
            }
        }
    }

    return encode_json({ error => 'Error: item not found' }) unless $id and $title;
    return encode_json(
        { timestamp => $timestamp, name => $title, id => $accn } );
}

sub upload_file {
    my %opts      = @_;
    my $filename  = '' . $FORM->param('input_upload_file');
    my $fh        = $FORM->upload('input_upload_file');

    #	print STDERR "upload_file: $filename\n";

    my $size = 0;
    my $path;
    if ($fh) {
        my $tmpfilename =
          $FORM->tmpFileName( $FORM->param('input_upload_file') );
        $path = 'upload/' . $filename;
        my $targetpath = $TEMPDIR . 'upload/';
        mkpath($targetpath);
        $targetpath .= $filename;

        #		print STDERR "temp files: $tmpfilename $targetpath\n";
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
	#print STDERR $user->user_name . ' ' . int($user->is_public) . "\n";
	return ($user && !$user->is_public);
}

sub load_genome {
    my %opts         = @_;
    my $name         = $opts{name};
    my $description  = $opts{description};
    my $link         = $opts{link};
    my $version      = $opts{version};
    my $type_id      = $opts{type_id};
    my $source_name  = $opts{source_name};
    my $restricted   = $opts{restricted};
    my $organism_id  = $opts{organism_id};
    my $user_name    = $opts{user_name};
    my $keep_headers = $opts{keep_headers};
    my $items        = $opts{items};

	#print STDERR "load_genome: organism_id=$organism_id name=$name description=$description version=$version type_id=$type_id restricted=$restricted\n";

	# Added EL: 7/8/2013.  Solves the problem when restricted is unchecked.
	# Otherwise, command-line call fails with '-organism_id' being passed to
	# restricted as option
    $restricted = ( $restricted && $restricted eq 'true' ) ? 1 : 0;

	# Check login
    if ( !$user_name || !$user->is_admin ) {
        $user_name = $user->user_name;
    }
    if ($user_name eq 'public') {
    	return encode_json({ error => 'Not logged in' });
    }

    # Check data items
    return encode_json({ error => 'No files specified' }) unless $items;
    $items = decode_json($items);
    my @files = map { catfile($TEMPDIR, $_->{path}) } @$items;

    # Setup staging area
    my $stagepath = catdir($TEMPDIR, 'staging');
    mkpath $stagepath;

    # Gather NCBI accession numbers if present
    my @accns;
    foreach my $item (@$items) {
        if ($item->{type} eq 'ncbi') {
            my $path = $item->{path};
            $path =~ s/\.\d+$//; # strip off version number
            push @accns, $path;
        }
    }

    # Submit workflow to add genome
    my ($workflow_id, $error_msg);
    if (@accns) { # NCBI accession numbers specified
        ($workflow_id, $error_msg) = create_genome_from_NCBI(
            user => $user,
            accns => \@accns
        );
    }
    else { # File-based load
        ($workflow_id, $error_msg) = create_genome_from_file(
            user => $user,
            metadata => {
                name => $name,
                description => $description,
                version => $version,
                source_name => $source_name,
                restricted => $restricted,
                organism_id => $organism_id,
                type_id => $type_id
            },
            files => \@files
        );
    }
    unless ($workflow_id) {
        return encode_json({ error => "Workflow submission failed: " . $error_msg });
    }

    # Get tiny link
    my $tiny_link = CoGe::Accessory::Web::get_tiny_link(
        url => $P->{SERVER} . "$PAGE_TITLE.pl?job_id=" . $workflow_id
    );

    return encode_json({ job_id => $workflow_id, link => $tiny_link });
}

sub get_load_log {
    my %opts         = @_;
    my $workflow_id = $opts{workflow_id};
    return unless $workflow_id;
    #TODO authenticate user access to workflow

    my (undef, $results_path) = get_workflow_paths($user->name, $workflow_id);
    return unless (-r $results_path);

    my $result_file = catfile($results_path, '1');
    return unless (-r $result_file);

    my $result = CoGe::Accessory::TDS::read($result_file);
    return unless $result;

    my $genome_id = (exists $result->{genome_id} ? $result->{genome_id} : undef);
    my $links = (exists $result->{links} ? $result->{links} : undef);

    return encode_json(
        {
            genome_id   => $genome_id,
            links       => $links
        }
    );
}

sub get_sequence_types {
    my $selected = 1;

    my $html;
    foreach my $type ( sort { $a->name cmp $b->name }
        $coge->resultset('GenomicSequenceType')->all() )
    {
        $html .= '<option value="' . $type->id . '"';
        if ( $selected && $type->id == $selected )
        {    #$type->name =~ /$selected/i) {
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

    my $type =
      $coge->resultset('GenomicSequenceType')
      ->find_or_create( { name => $name, description => $desc } );
    return unless $type;

    return $type->id;
}

sub create_organism {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    return unless $name and $desc;
    #print STDERR "create_organism $name $desc\n";

    my $organism =
      $coge->resultset('Organism')
      ->find_or_create( { name => $name, description => $desc } );
    return unless $organism;

    return $organism->id;
}

sub search_organisms {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @organisms = $coge->resultset("Organism")->search(
        \[
            'name LIKE ? OR description LIKE ?'
            ,    #FIXME security hole: need to check 'restricted'
            [ 'name', $search_term ], [ 'description', $search_term ]
        ]
    );

    # Limit number of results displayed
    if ( @organisms > $MAX_SEARCH_RESULTS ) {
        return encode_json( { timestamp => $timestamp, items => undef } );
    }

    my @results;
    foreach ( sort { $a->name cmp $b->name } @organisms ) {
        push @results, { 'label' => $_->name, 'value' => $_->id };
    }

    return encode_json( { timestamp => $timestamp, items => \@results } );
}

sub search_users {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @users = $coge->resultset("User")->search(
        \[
            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
            [ 'user_name',  $search_term ],
            [ 'first_name', $search_term ],
            [ 'last_name',  $search_term ]
        ]
    );

    # Limit number of results displayed
    # if (@users > $MAX_SEARCH_RESULTS) {
    # 	return encode_json({timestamp => $timestamp, items => undef});
    # }

    return encode_json(
        {
            timestamp => $timestamp,
            items     => [ sort map { $_->user_name } @users ]
        }
    );
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

sub send_error_report {
    my %opts = @_;
    my $load_id = $opts{load_id};
    my $job_id = $opts{job_id};

    my @paths= ($P->{SECTEMPDIR}, $PAGE_TITLE, $user->name, $load_id, "staging");

    # Get the staging directory
    my $staging_dir = File::Spec->catdir(@paths);

    my $url = $P->{SERVER} . "$PAGE_TITLE.pl?";
    $url .= "job_id=$job_id;" if $job_id;
    $url .= "load_id=$load_id";

    my $email = $P->{SUPPORT_EMAIL};

    my $body =
        "Load failed\n\n"
        . 'For user: '
        . $user->name . ' id='
        . $user->id . ' '
        . $user->date . "\n\n"
        . "staging_directory: $staging_dir\n\n"
        . "tiny link: $url\n\n";
    #$body .= get_load_log();

    CoGe::Accessory::Web::send_email(
        from    => $email,
        to      => $email,
        subject => "Load error notification from $PAGE_TITLE",
        body    => $body
    );
}
