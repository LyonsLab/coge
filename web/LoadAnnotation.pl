#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadGenome.pl, replicate changes when applicable.

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path);
use HTML::Template;
use JSON::XS;
use URI::Escape::JavaScript qw(escape unescape);
use File::Path;
use File::Copy;
use File::Basename;
use File::Slurp;
use File::Spec::Functions qw(catdir catfile);
use File::Touch;
use File::Listing qw(parse_dir);
use Sort::Versions;
no warnings 'redefine';

use vars qw(
  $CONF $PAGE_TITLE $LINK $EMBED $TEMPDIR $USER $DB $FORM
  %FUNCTION $LOAD_ID $WORKFLOW_ID
);

$PAGE_TITLE = 'LoadAnnotation';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

# Get workflow_id and load_id for previous load if specified.  Otherwise
# generate a new load_id for data upload.
$WORKFLOW_ID = $FORM->Vars->{'wid'} || $FORM->Vars->{'job_id'}; # wid is new name, job_id is legacy name
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = get_upload_path($USER->name, $LOAD_ID);

$EMBED = $FORM->param('embed');

%FUNCTION = (
    upload_file     => \&upload_file,
    get_sources     => \&get_sources,
    create_source   => \&create_source,
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
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param( PAGE_TITLE => $PAGE_TITLE,
					      TITLE      => "LoadAnnotation",
        				  PAGE_LINK  => $LINK,
        				  HOME       => $CONF->{SERVER},
                          HELP       => 'LoadAnnotation',
                          WIKI_URL   => $CONF->{WIKI_URL} || '',
                          USER       => $USER->display_name || '' );
        $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
        
        my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
        $link = CoGe::Accessory::Web::get_tiny_link( url => $link );
    
        $template->param( ADMIN_ONLY => $USER->is_admin,
                          CAS_URL    => $CONF->{CAS_URL} || '',
                          COOKIE_NAME => $CONF->{COOKIE_NAME} || '' );
    }
    
    $template->param( BODY => generate_body() );
    return $template->output;
}

sub generate_body {
    if ( $USER->user_name eq 'public' ) {
        my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
        $template->param( LOGIN     => 1 );
        return $template->output;
    }

    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );

    my $gid;
    $gid = $FORM->param('gid') if defined $FORM->param('gid');
    if ($gid) {
        my $genome = $DB->resultset('Genome')->find($gid);
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
        API_BASE_URL  => $CONF->{SERVER} . 'api/v1/', #TODO move into config file or module
        HELP_URL      => 'https://genomevolution.org/wiki/index.php/LoadAnnotation',
        SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
        DEFAULT_TAB              => 0,
        FILE_SELECT_SINGLE       => 1,
        DISABLE_IRODS_GET_ALL    => 1,
        DISABLE_FTP_GET_ALL      => 1,
        MAX_IRODS_LIST_FILES     => 1000,
        MAX_IRODS_TRANSFER_FILES => 30,
        MAX_FTP_FILES            => 30,
        USER                     => $USER->user_name,
    );
    $template->param( SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
                      SPLASH_CONTENTS    => 'This page allows you to load genome annotation in GFF file format.' );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

sub upload_file {
    my %opts      = @_;
    my $filename  = '' . $FORM->param('input_upload_file');
    my $fh        = $FORM->upload('input_upload_file');
    #print STDERR "upload_file: $filename\n";

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
        touch($targetpath . '.done'); # for JEX
        $size = -s $fh;

        # mdb added 2/9/17 -- delete temp file
        unlink($tmpfilename) or warn "Couldn't delete file '$tmpfilename': $!";
    }

    return encode_json(
        {
            filename  => $filename,
            path      => $path,
            size      => $size
        }
    );
}

sub get_sources {
    #my %opts = @_;

    my %unique;
    foreach ( $DB->resultset('DataSource')->all() ) {
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
      $DB->resultset('DataSource')
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

    my $url = $CONF->{SERVER} . "$PAGE_TITLE.pl?";
    $url .= "job_id=$job_id;" if $job_id;
    $url .= "load_id=$load_id";

    my $email = $CONF->{SUPPORT_EMAIL};

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
