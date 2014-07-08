#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use Digest::MD5 qw(md5_base64);
use HTML::Template;
use Data::Dumper;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $COOKIE_NAME $coge);
$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};

# set this to 1 to print verbose messages to logs
$DEBUG   = 0;
$TEMPDIR = $P->{TEMPDIR};
$TEMPURL = $P->{TEMPURL};
$|       = 1;               # turn off buffering
$DATE    = sprintf(
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

my $pj = new CGI::Ajax( gen_html => \&gen_html, );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
#print $pj->build_html($FORM, \&gen_html);
print "Content-Type: text/html\n\n" . gen_html();

sub gen_html {
    my $html;    # =  "Content-Type: text/html\n\n";
    my ( $body, $seq_names, $seqs ) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'MSAView.tmpl' );

    $template->param( TITLE => 'Multiple Sequence Alignment Viewer' );
    #    $template->param(PAGE_TITLE=>'MSAView');
    $template->param( HELP => "/wiki/index.php?title=MSAView" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE => $DATE );
    $template->param( LOGO_PNG  => "MSAView-logo.png" );
    $template->param( BODY      => $body );
    $template->param( SEQ_NAMES => $seq_names ) if $seq_names;
    $template->param( SEQS      => $seqs ) if $seqs;
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $html;
    my ( $seq_names, $seqs );
    if ( $FORM->param('file') ) {
        my $fh = $FORM->upload('file');
        my $content;
        while (<$fh>) { $content .= $_ }
        ( $html, $seq_names, $seqs ) = gen_alignment_view($content);
    }
    else {
        $html .= qq{
<form method="post" enctype="multipart/form-data">
<DIV>Please select an alignment file in fasta format to upload:
<input type = "file" name= "file">
</DIV>
<input type = "submit" name="GO" value="GO">
</FORM>
}
    }
    return ( $html, $seq_names, $seqs );
}

sub gen_alignment_view {
    my $content = shift;
    my ( $seqs, $order, $stats ) = parse_fasta($content);
    my $html;
    my ( $seq_names, $out_seqs );
    $seq_names = join "\n", @$order;
    $out_seqs = join "\n", map { $seqs->{$_} } @$order;
    $html .= "<hr>";
    $html .= "Stats:";
    $html .= "<table>";
    $html .= "<tr>";
    $html .= "<td><div class=small>";
    $html .= "Number of sequences";
    $html .= "</div>";
    $html .= "<td><div class=small>";
    $html .= $stats->{num_seqs};
    $html .= "</div>";

    $html .= "<tr>";
    $html .= "<td><div class=small>";
    $html .= "Alignment length";
    $html .= "</div>";
    $html .= "<td><div class=small>";
    $html .= $stats->{length} . " characters";
    $html .= "</div>";
    foreach my $num (qw(100 75 50 25 5)) {
        $html .= "<tr>";
        $html .= "<td><div class=small>";
        $html .= "Characters with " . $num . "% identity";
        $html .= "</div>";
        $html .= "<td><div class=small>";
        my $stat = $stats->{percent}{$num} || 0;
        my $val = sprintf( "%.2f", $stat / $stats->{length} ) * 100;
        $html .= $stat . "(" . $val . "%)";
        $html .= "</div>";
    }
    $html .= "</table>";
    return ( $html, $seq_names, $out_seqs );
}

sub parse_fasta {
    my $text = shift;
    my %seqs;
    my @order;
    foreach my $ent ( split /\n>/, $text ) {
        $ent =~ s/>//g;
        my ( $name, $seq ) = split /\n/, $ent, 2;
        $seq =~ s/\n//g;
        $seqs{$name} .= $seq;
        push @order, $name;
    }
    my ( $cons, $stats ) = generate_consensus_seq( [ values %seqs ] );
    $seqs{"consensus"} = $cons;
    push @order, "consensus";
    return ( \%seqs, \@order, $stats );
}

sub generate_consensus_seq {
    my $seqs     = shift;
    my $num_seqs = scalar @$seqs;
    return unless $num_seqs;
    my $seq_len = length $seqs->[0];
    my $con;
    my %stats;
    $stats{length}   = $seq_len;
    $stats{num_seqs} = $num_seqs;
    for ( my $i = 0 ; $i < $seq_len ; $i++ ) {
        my $let;
        my %res;
        for ( my $j = 0 ; $j < $num_seqs ; $j++ ) {
            my $chr = substr $seqs->[$j], $i, 1;
            $res{$chr}++ unless $chr eq "-";
            $let = $chr unless $let;
            $let = " " unless $chr eq $let;
        }

        my ($max_chr) = sort { $b <=> $a } values %res;
        $stats{percent}{100}++ unless $let eq " ";
        my $pid = $max_chr / $num_seqs;
        $stats{percent}{75}++ if $pid >= 0.75;
        $stats{percent}{50}++ if $pid >= 0.5;
        $stats{percent}{25}++ if $pid >= 0.25;
        $stats{percent}{5}++  if $pid >= 0.05;
        $let = "*" if $pid >= .75 && $let eq " ";
        $let = ":" if $pid >= .5  && $let eq " ";
        $let = "." if $pid >= .25 && $let eq " ";
        $con .= $let;
    }
    return $con, \%stats;
}
