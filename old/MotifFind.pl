#! /usr/bin/perl -w

use strict;
use CGI;

#use CGI::Ajax;
use JSON::XS;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use URI::Escape;
use Data::Dumper;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL);
$P = CoGe::Accessory::Web::get_defaults();

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM = new CGI;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL         = $P->{URL};
$COGEDIR     = $P->{COGEDIR};
$TEMPDIR     = $P->{TEMPDIR} . "GenomeList/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "GenomeList/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    ticket   => $cas_ticket,
    coge     => $coge,
    this_url => $FORM->url()
) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

%FUNCTION = (
    gen_html   => \&gen_html,
    gen_go_run => \&gen_go_run,
);
dispatch();

sub dispatch {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};
    if ($fname) {

        #my %args = $FORM->Vars;
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
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => '/wiki/index.php?title=BLANK' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER       => $name );
    $template->param( TITLE      => qq{Motif Find} );
    $template->param( PAGE_TITLE => qq{Motif Find} );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE       => $DATE );
    $template->param( BODY       => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_params {
    my $params;
    $params .= qq{
'args__motif_choice','args__'+pageObj.motif_choice,
'args__gene_names', 'gene_names',
'args__search_seqs', 'search_seqs',
'args__motif_names', 'motif_names',
'args__up_stream', 'up_stream',
'args__down_stream', 'down_stream',

};
    return $params;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'MotifFind.tmpl' );
    $template->param( PAGE_NAME => $FORM->url );

    $template->param( MOTIF_SELECT => motif_list() );
    return $template->output;
}

sub gen_go_run {
    my %opts = &gen_params();

    print STDERR %opts;

    #my $motif_choice = $opts{motif_choice};

    #my $output = "Test worked!\n";
    #$output .= "<pre>"."Args: ".Dumper (\%opts)."</pre>";
}

sub motif_list {
    my $MOTIFFILE = "/opt/apache/CoGe-Dev/bin/MotifView/newMotifHash";
    my $motifhash = do $MOTIFFILE || print STDERR "Cannot open $MOTIFFILE";
    my @opts;

    my $html;

    for my $mot (
        sort { $motifhash->{$a}{'consensus'} cmp $motifhash->{$b}{'consensus'} }
        keys %$motifhash
      )
    {
        my $seq            = $motifhash->{$mot}{'consensus'};
        my $motifrealcolor = $motifhash->{$mot}{'color'};
        my $motifrgbcolor  = $motifhash->{$mot}{'rgb'};
        my $motifhexcolor  = $motifhash->{$mot}{'hex'};
        my $auth           = $motifhash->{$mot}{'auth'};
        my $litref         = $motifhash->{$mot}{'litRef'};
        my $title          = $motifhash->{$mot}{'title'};
        my $motifname      = $motifhash->{$mot}{'name'};
        my $iupac          = $motifhash->{$mot}{'iupac'};
        my $family         = $motifhash->{$mot}{'familyname'};
        my $source         = $motifhash->{$mot}{'source'};
        my $stress         = $motifhash->{$mot}{'stress'};
        my $hex            = $motifhash->{$mot}{'hex'};

        $seq            =~ s/\"//g;
        $motifrealcolor =~ s/\"//g;
        $motifrgbcolor  =~ s/\"//g;
        $motifhexcolor  =~ s/\"//g;
        $auth           =~ s/\"//g;
        $litref         =~ s/\"//g;
        $title          =~ s/\"//g;

        $motifname =~ s/\"//g;
        $motifname =~ s/\,/_/g;
        $motifname =~ s/\-/_/g;
        $motifname =~ s/\s+/_/g;

        $motifrealcolor =~ s/\s+//g;
        $motifrgbcolor  =~ s/\s+//g;
        $motifhexcolor  =~ s/\s+//g;
        $auth           =~ s/\"//g;
        $auth           =~ s/:/ /g;
        $auth           =~ s/\s+/_/g;
        $litref         =~ s/\"//g;
        $litref         =~ s/:/ /g;
        $litref         =~ s/\s+/_/g;
        $title          =~ s/\"//g;
        $title          =~ s/:/ /g;
        $title          =~ s/\s+/_/g;
        $seq            =~ s/\s+//g;
        $title          =~ s/\-//g;
        chomp $title;

        my $showmotifname = $motifname;
        if ( $motifname =~ /\;/ ) {
            my $motifnamelength = length($showmotifname);
            $showmotifname = substr $motifname, 0, 10 if $motifnamelength > 11;
        }

        my $motifval =
"$mot\|$family\|$motifname\|$seq\|$iupac\|$auth\|$title\|$litref\|NA\|\|NA\|$stress\|$motifrealcolor\|$motifrgbcolor\|$hex\|$source\|NA\|NA";

        if ( $seq =~ /\w+/ ) {
            push @opts,
                "<OPTION value=\"$motifval\">"
              . $iupac . " : "
              . $showmotifname
              . "</OPTION>";
        }
    }
    $html .=
        qq{<FONT CLASS ="small" id="motif_count">Motif count: }
      . scalar @opts
      . qq{</FONT>\n<BR>\n};
    $html .=
qq{<SELECT class="backbox" id="motif_select" SIZE="15" MULTIPLE onclick="show_add();" ondblclick="add_selected_motifs();">\n};

    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";

    $html =~ s/OPTION/OPTION SELECTED/;

    return $html;
}
