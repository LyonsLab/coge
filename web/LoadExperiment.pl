#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadGenome.pl, replicate changes when applicable.

use strict;

use CGI;
use HTML::Template;
use JSON::XS;
use URI::Escape::JavaScript qw(unescape);
use File::Path;
use File::Copy;
use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use File::Listing qw(parse_dir);
use File::Slurp;
use LWP::Simple;
use URI;
use Sort::Versions;
use Data::Dumper;

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::IRODS;
use CoGe::Accessory::TDS;
use CoGe::Accessory::Utils;
use CoGe::Core::Genome qw(genomecmp);
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path get_irods_file);

no warnings 'redefine';

use vars qw(
  $P $PAGE_TITLE $TEMPDIR $USER $coge $FORM $LINK $EMBED
  %FUNCTION $MAX_SEARCH_RESULTS $CONFIGFILE $LOAD_ID $WORKFLOW_ID
);

$PAGE_TITLE = 'LoadExperiment';

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$CONFIGFILE = catfile($ENV{COGE_HOME}, 'coge.conf');

# Get workflow_id and load_id for previous load if specified.  Otherwise
# generate a new load_id for data upload.
$WORKFLOW_ID = $FORM->Vars->{'wid'} || $FORM->Vars->{'job_id'}; # wid is new name, job_id is legacy name
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = get_upload_path($USER->name, $LOAD_ID);

$EMBED = $FORM->param('embed');

$MAX_SEARCH_RESULTS = 1000;

%FUNCTION = (
    irods_get_path          => \&irods_get_path,
    irods_get_file          => \&irods_get_file,
    load_from_ftp           => \&load_from_ftp,
    ftp_get_file            => \&ftp_get_file,
    upload_file             => \&upload_file,
    get_sources             => \&get_sources,
    create_source           => \&create_source,
    search_genomes          => \&search_genomes,
    search_users            => \&search_users,
    send_error_report       => \&send_error_report
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub generate_html {
    # Check for finished result # mdb removed 3/4/15 no longer auto-redirect, make user select result
#    if ($LOAD_ID) {
#        my $log = get_load_log(workflow_id => $WORKFLOW_ID);
#        if ($log) {
#            my $res = decode_json($log);
#            if ($res->{experiment_id}) {
#                my $url = 'ExperimentView.pl?eid=' . $res->{experiment_id};
#                print $FORM->redirect(-url => $url);
#            }
#        }
#    }
    
    my $template;
    
    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param( PAGE_TITLE => $PAGE_TITLE,
		                  TITLE      => "Load Experiment",
        	              PAGE_LINK  => $LINK,
			              HOME       => $P->{SERVER},
                          HELP       => 'LoadExperiment',
                          WIKI_URL   => $P->{WIKI_URL} || '',
			              ADJUST_BOX => 1,
                          ADMIN_ONLY => $USER->is_admin,
                          USER       => $USER->display_name || '',
                          CAS_URL    => $P->{CAS_URL} || ''
        );
        $template->param( LOGON      => 1 ) unless $USER->is_public;
    }

    $template->param( BODY => generate_body() );
    return $template->output;
}

sub generate_body {
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
    
    # Force login
    if ( $USER->is_public ) {
        $template->param( LOGIN => 1 );
        return $template->output;
    }

    # Set genome ID if specified
    my $gid = $FORM->param('gid');
    if ($gid) {
        my $genome = $coge->resultset('Genome')->find($gid);
        if ($genome && $USER->has_access_to_genome($genome)) { # check permission
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
        HELP_URL      => 'https://genomevolution.org/wiki/index.php/LoadExperiment',
        SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
        DEFAULT_TAB              => 0,
        MAX_IRODS_LIST_FILES     => 1000,
        MAX_IRODS_TRANSFER_FILES => 30,
        MAX_FTP_FILES            => 30,
        USER                     => $USER->user_name
    );
    $template->param( SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
                      SPLASH_CONTENTS    => 'This page allows you to load quantitative, polymorphism, or alignment data onto a genome from a variety of file formats.' );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

#sub irods_get_path {
#    my %opts      = @_;
#    my $path      = $opts{path};
#    $path = unescape($path);
#    #print STDERR "irods_get_path: $path\n";
#    my $username = $USER->name;
#    my $basepath = $P->{IRODSDIR};
#    $basepath =~ s/\<USER\>/$username/;
#    $path = $basepath unless $path;
#
#    if ( $path !~ /^$basepath/ ) {
#        print STDERR "Attempt to access '$path' denied (basepath='$basepath')\n";
#        return;
#    }
#
#    my $result = CoGe::Accessory::IRODS::irods_ils($path, escape_output => 1);
#    #print STDERR "irods_get_path ", Dumper $result, "\n";
#    my $error  = $result->{error};
#    if ($error) {
#        my $email = $P->{SUPPORT_EMAIL};
#        my $body =
#            "irods ils command failed\n\n"
#          . 'User: '
#          . $USER->name . ' id='
#          . $USER->id . ' '
#          . $USER->date . "\n\n"
#          . $error . "\n\n"
#          . $P->{SERVER};
#        CoGe::Accessory::Web::send_email(
#            from    => $email,
#            to      => $email,
#            subject => "System error notification from $PAGE_TITLE",
#            body    => $body
#        );
#        return encode_json( { error => $error } );
#    }
#    return encode_json( { path => $path, items => $result->{items} } );
#}
#
#sub irods_get_file {
#    my %opts = @_;
#    my $path = $opts{path};
#    
#    $path = unescape($path);
#    
#    my $result = get_irods_file($path, $TEMPDIR);
#
#    return encode_json( { path => $result->{path}, size => $result->{size} } );
#}

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

    my $path         = catdir('ftp', $filepath, $filename);
    my $fullfilepath = catdir($TEMPDIR, 'ftp', $filepath);
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

    #$request->content_type("text/xml; charset=utf-8"); # mdb removed 3/30/15 COGE-599
    my $response = $ua->request($request);
    #print STDERR "content: <begin>", $response->content , "<end>\n"; # debug
    if ( $response->is_success() ) {
        #my $header = $response->header;
        my $result = $response->content;
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
                error     => "Failed: $status"
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

sub ncbi_search {
    my %opts      = @_;
    my $accn      = $opts{accn};
    my $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$accn";
    my $result = get($esearch);

    #print STDERR $result;

    my $record = XMLin($result);

    #print STDERR Dumper $record;

    my $id = $record->{IdList}->{Id};
    print STDERR "id = $id\n";

    my $title;
    if ($id) {
        $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=$id";
        my $result = get($esearch);
        #print STDERR $result;
        $record = XMLin($result);
        #print STDERR Dumper $record;

        foreach ( @{ $record->{DocSum}->{Item} } )
        {    #FIXME use grep here instead
            if ( $_->{Name} eq 'Title' ) {
                $title = $_->{content};
                print STDERR "title=$title\n";
                last;
            }
        }
    }

    return unless $id and $title;
    return encode_json(
        { name => $title, id => $id } );
}

sub upload_file {
    my %opts      = @_;
    my $filename  = '' . $FORM->param('input_upload_file');
    my $fh        = $FORM->upload('input_upload_file');

    #	print STDERR "upload_file: $filename\n";

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

#sub search_genomes
#{    # FIXME: common with LoadAnnotation et al., move into web service
#    my %opts        = @_;
#    my $search_term = $opts{search_term};
#    my $timestamp   = $opts{timestamp};
#    #print STDERR "$search_term $timestamp\n";
#    return unless $search_term;
#
#    # Perform search
#    my $id = $search_term;
#    $search_term = '%' . $search_term . '%';
#
#    # Get all matching organisms
#    my @organisms = $coge->resultset("Organism")->search(
#        \[
#            'name LIKE ? OR description LIKE ?',
#            [ 'name',        $search_term ],
#            [ 'description', $search_term ]
#        ]
#    );
#
#    # Get all matching genomes
#    my @genomes = $coge->resultset("Genome")->search(
#        \[
#            'genome_id = ? OR name LIKE ? OR description LIKE ?',
#            [ 'genome_id',   $id ],
#            [ 'name',        $search_term ],
#            [ 'description', $search_term ]
#        ]
#    );
#
#    # Combine matching genomes with matching organism genomes, preventing duplicates
#    my %unique;
#    map {
#        $unique{ $_->id } = $_ if ( $USER->has_access_to_genome($_) )
#    } @genomes;
#    foreach my $organism (@organisms) {
#        map {
#            $unique{ $_->id } = $_ if ( $USER->has_access_to_genome($_) )
#        } $organism->genomes;
#    }
#
#    # Limit number of results displayed
#    if ( keys %unique > $MAX_SEARCH_RESULTS ) {
#        return encode_json( { timestamp => $timestamp, items => undef } );
#    }
#
#    my @items;
#    #print STDERR Dumper \@items, "\n";
#    foreach ( sort genomecmp values %unique ) {    #(keys %unique) {
#        push @items, { label => $_->info, value => $_->id };
#    }
#
#    return encode_json( { timestamp => $timestamp, items => \@items } );
#}

#sub search_users {
#    my %opts        = @_;
#    my $search_term = $opts{search_term};
#    my $timestamp   = $opts{timestamp};
#
#    #print STDERR "$search_term $timestamp\n";
#    return unless $search_term;
#
#    # Perform search
#    $search_term = '%' . $search_term . '%';
#    my @users = $coge->resultset("User")->search(
#        \[
#            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
#            [ 'user_name',  $search_term ],
#            [ 'first_name', $search_term ],
#            [ 'last_name',  $search_term ]
#        ]
#    );
#
#    # Limit number of results displayed
#    # if (@users > $MAX_SEARCH_RESULTS) {
#    # 	return encode_json({timestamp => $timestamp, items => undef});
#    # }
#
#    return encode_json(
#        {
#            timestamp => $timestamp,
#            items     => [ sort map { $_->user_name } @users ]
#        }
#    );
#}

sub get_sources {
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
    unless ($load_id and $job_id) {
        print STDERR "LoadExperiment::send_error_report: missing required params\n";
        return;
    }

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

    my $log = get_debug_log(workflow_id => $job_id);
    $body .= $log if $log;

    CoGe::Accessory::Web::send_email(
        from    => $email,
        to      => $email,
        subject => "Load error notification from $PAGE_TITLE",
        body    => $body
    );
}
