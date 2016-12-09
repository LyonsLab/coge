#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadExperiment.pl & LoadAnnotation.pl, replicate changes when applicable.

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage qw(get_workflow_paths get_upload_path);
use HTML::Template;
use LWP::Simple qw(get);
use XML::Simple qw(XMLin);
use JSON::XS;
use Sort::Versions;
use File::Path qw(mkpath);
use File::Copy qw(copy);
use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use File::Touch;
use File::Listing qw(parse_dir);
use File::Slurp;
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;
no warnings 'redefine';

use vars qw(
  $CONF $PAGE_TITLE $TEMPDIR $USER $DB $FORM $LINK $WORKFLOW_ID $EMBED
  %FUNCTION $LOAD_ID
);

$PAGE_TITLE = 'LoadGenome';

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
    upload_file          => \&upload_file,
    search_ncbi_nucleotide => \&search_ncbi_nucleotide,
    get_sequence_types   => \&get_sequence_types,
    create_sequence_type => \&create_sequence_type,
    create_source        => \&create_source,
    get_sources          => \&get_sources,
	send_error_report    => \&send_error_report
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub generate_html {
    # Check for finished result # mdb removed 3/4/15 no longer auto-redirect, make user select result
#    if ($WORKFLOW_ID) {
#    	my $log = get_load_log(workflow_id => $WORKFLOW_ID);
#    	if ($log) {
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
					      TITLE      => "Load Genome",
        				  PAGE_LINK  => $LINK,
					      HOME       => $CONF->{SERVER},
                          HELP       => 'LoadGenome',
                          WIKI_URL   => $CONF->{WIKI_URL} || '',
                          USER       => $USER->display_name || ''
        );
        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
        my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
        $link = CoGe::Accessory::Web::get_tiny_link( url => $link );
    
        $template->param( ADMIN_ONLY  => $USER->is_admin,
                          CAS_URL     => $CONF->{CAS_URL} || '',
                          COOKIE_NAME => $CONF->{COOKIE_NAME} || ''
        );
    }
    
    $template->param( BODY => generate_body() );
    return $template->output;
}

sub generate_body {
    if ( $USER->user_name eq 'public' ) {
        my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param(
            PAGE_NAME => "$PAGE_TITLE.pl",
            LOGIN     => 1
        );
        return $template->output;
    }
    
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        MAIN          => 1,
        PAGE_TITLE    => $PAGE_TITLE,
        EMBED         => $EMBED,
        PAGE_NAME     => $PAGE_TITLE . '.pl',
        LOAD_ID       => $LOAD_ID,
        WORKFLOW_ID   => $WORKFLOW_ID,
        API_BASE_URL  => $CONF->{SERVER} . 'api/v1/', #TODO move into config file or module
        HELP_URL      => 'https://genomevolution.org/wiki/index.php/LoadGenome',
        SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
        ENABLE_NCBI              => 1,
        DEFAULT_TAB              => 0,
        MAX_IRODS_LIST_FILES     => 1000,
        MAX_IRODS_TRANSFER_FILES => 30,
        MAX_FTP_FILES            => 50,
        USER                     => $USER->user_name
    );
    
    $template->param( SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
                      SPLASH_CONTENTS    => 'This page allows you to load genome sequences in FASTA file format.' );

    my $oid = $FORM->param("oid");
    my $organism = $DB->resultset('Organism')->find($oid) if $oid;
    if ($organism) {
        $template->param(ORGANISM_NAME => $organism->name);
    }

    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;

    return $template->output;
}

sub search_ncbi_nucleotide { #TODO this can be done client-side instead, see Entrez.js
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
    my $upload_file = $FORM->param('input_upload_file');
    my $filename  = '' . $upload_file;
    my $fh        = $FORM->upload('input_upload_file');

    #   print STDERR "upload_file: $filename\n";

    my $size = 0;
    my $path;
    if ($fh) {
        my $tmpfilename = $FORM->tmpFileName( $upload_file );
        $path = catfile('upload', $filename);
        my $targetpath = catdir($TEMPDIR, 'upload');
        mkpath($targetpath);
        $targetpath = catfile($targetpath, $filename);

        #print STDERR "temp files: $tmpfilename $targetpath\n";
        copy( $tmpfilename, $targetpath );
        touch($targetpath . '.done'); # for JEX
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

sub get_sequence_types {
    my $selected = 1;

    my $html;
    foreach my $type ( sort { $a->name cmp $b->name }
        $DB->resultset('GenomicSequenceType')->all() )
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
      $DB->resultset('GenomicSequenceType')
      ->find_or_create( { name => $name, description => $desc } );
    return unless $type;

    return $type->id;
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

    my $source = $DB->resultset('DataSource')->find_or_create(
        { name => $name, description => $desc, link => $link } );
    return unless ($source);

    return $name;
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
