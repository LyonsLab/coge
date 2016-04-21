#! /usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use CGI;
use DBI;
use Data::Dumper;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use File::Path;
use Digest::MD5 qw(md5_base64);
use Benchmark;
use DBIxProfiler;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM $URL $TEMPURL $MAUVE $COGE_MAUVE $MAUVE_MATRIX $COOKIE_NAME);
$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE      = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);
$PAGE_NAME = "GenomeAlign.pl";

$TEMPDIR = $P->{TEMPDIR} . "GenomeAlign/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL      = $P->{TEMPURL} . "GenomeAlign/";
$FORM         = new CGI;
$MAUVE        = $P->{MAUVE};
$COGE_MAUVE   = $P->{COGE_MAUVE};
$MAUVE_MATRIX = $P->{MAUVE_MATRIX};
$DBNAME       = $P->{DBNAME};
$DBHOST       = $P->{DBHOST};
$DBPORT       = $P->{DBPORT};
$DBUSER       = $P->{DBUSER};
$DBPASS       = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    cookie_name => $COOKIE_NAME,
    ticket      => $cas_ticket,
    coge        => $coge,
    this_url    => $FORM->url()
) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

$SIG{'__WARN__'} = sub { };    #silence warnings

my %ajax = CoGe::Accessory::Web::ajax_func();

my %FUNCTION = (
    gen_html      => \&gen_html,
    gen_data      => \&gen_data,
    run_alignment => \&run_alignment,
    %ajax,
);

#my $pj = new CGI::Ajax(%FUNCTION);
#$pj->js_encode_function('escape');
#my $t1 = new Benchmark;
if ( $FORM->param('jquery_ajax') ) {
    dispatch();
}
else {
    print $FORM->header, "\n", gen_html();
}

sub dispatch {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};
    if ($fname) {
        if ( $args{args} ) {
            my @args_list = split( /,/, $args{args} );
            print $FORM->header, $FUNCTION{$fname}->(@args_list);
        }
        else {
            print $FORM->header, $FUNCTION{$fname}->(%args);
        }
    }
}

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => 'GenomeAlign' );
    $template->param( HELP       => '/wiki/index.php?title=GenomeAlign' );

    # print STDERR "user is: ",$USER,"\n";
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER     => $name );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE     => $DATE );
    my $list_name = $FORM->param('list_name') || $FORM->param('ln');
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GenomeAlign.tmpl' );
    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $sort_by_type     = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $prefs =
      CoGe::Accessory::Web::load_settings( user => $USER, page => $PAGE_NAME );
    $prefs = {} unless $prefs;
    $template->param( 'TEMPURL' => $TEMPURL );
    my $dsgids = [];
    $dsgids = read_file() if $BASEFILE;    #: $opts{feature_list};

    foreach my $item ( $form->param('dsgid') ) {
        foreach my $item2 ( split /(::)|(,)/, $item ) {
            push @$dsgids, $item2 if $item2 =~ /^\d+_?\d*$/;
        }
    }
    my ( $table, $count ) = generate_table( dsgids => $dsgids );
    $template->param( 'GENOME_COUNT' => $count );
    my $genomelist_link = "GenomeList.pl?dsgid=" . join( ",", @$dsgids );
    $template->param( 'GENOMELIST_LINK' => $genomelist_link );
    if ($table) {
        $template->param( INFO => $table );
        return $template->output;
    }
    else {
        return "No genomes were specified.";
    }
}

sub generate_table {
    my %opts   = @_;
    my $dsgids = $opts{dsgids};
    return unless @$dsgids;
    my @table;
    my $count = 1;
    foreach my $dsgid (@$dsgids) {
        my $dsg  = $coge->resultset('Genome')->find($dsgid);
        my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
        my $desc = $dsg->description ? $dsg->description : join(
            "; ",
            map {
                qq{<span class=link onclick=window.open('OrganismView.pl?org_desc=$_')>$_</span>}
              } split /;\s*/,
            $dsg->organism->description
        );

    	my $chr_count = $dsg->chromosome_count;
        my $length    = $dsg->length;
        my $type      = $dsg->type->name;
        push @table, {
            COUNT     => $count,
            DSGID     => $dsgid,
            NAME      => $name,
            DESC      => $desc,
            VER       => $dsg->version,
            TYPE      => $type,
            CHR_COUNT => commify($chr_count),
            LENGTH    => commify($length),

        };
        $count++;
    }
    $count--;
    return \@table, $count;
}

sub run_alignment {
    my %opts     = @_;
    my $dsgids   = $opts{dsgids};
    my $basename = $opts{basename};
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    $dsgids =~ s/^,+//;
    $dsgids =~ s/,+$//;
    my $cmd =
      $COGE_MAUVE . " -mauve_bin $MAUVE -matrix $MAUVE_MATRIX -dsgid '$dsgids'";
    my $outfile = $cogeweb->basefile . ".mauve.aln";
    $cmd .= " -out_file " . $outfile;
    CoGe::Accessory::Web::write_log( "Running $cmd", $cogeweb->logfile );
    my $output;
    open( CMD, "$cmd |" );
    open( OUT, ">" . $cogeweb->logfile );

    while (<CMD>) {
        print OUT $_;
    }
    close CMD;
    close OUT;
    CoGe::Accessory::Web::write_log( "########", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "LOG FROM $COGE_MAUVE",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "########", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "$output",  $cogeweb->logfile );

    my $logfile = $cogeweb->logfile;
    $outfile =~ s/$TEMPDIR/$TEMPURL/;
    $logfile =~ s/$TEMPDIR/$TEMPURL/;
    my $html;
    $html .= qq{
<div>Results:</div>
<div><a target=_new href=$outfile>Alignment File</a></div>
<div><a target=_new href=$logfile>Log File</a></div>
};
    return ($html);
}

sub gen_data {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
}
