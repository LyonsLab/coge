#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGeX;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use DBI;
no warnings 'redefine';

# for security purposes
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $COOKIE_NAME $coge);
$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$TEMPDIR   = $P->{TEMPDIR};
$TEMPURL   = $P->{TEMPURL};

# set this to 1 to print verbose messages to logs
$DEBUG = 0;

$|    = 1;         # turn off buffering
$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM                 = new CGI;
$CGI::POST_MAX        = 60 * 1024 * 1024;    # 24MB
$CGI::DISABLE_UPLOADS = 0;

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
($USER) = CoGe::Accessory::Web->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

my $pj = new CGI::Ajax( gen_html => \&gen_html, );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');

#print $pj->build_html($FORM, \&gen_html);
print "Content-Type: text/html\n\n";
print gen_html();

sub gen_html {
    my $form = shift || $FORM;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'HSP Viewer' );
    $template->param( PAGE_TITLE => 'HSPView' );
    #$template->param( HELP       => "/wiki/index.php?title=HSPView" );
    $template->param( HELP       => $P->{SERVER} );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE => $DATE );
    my $report_file = $form->param('blast_report') || $form->param('report');
    $report_file = $TEMPDIR . "/" . $report_file
      unless $report_file =~ /$TEMPDIR/;
    my $dbfile = $form->param('db');
    if ( $report_file && -r $report_file ) {
        my $hsp_num = $form->param('hsp_num') || $form->param('num');
        $report_file =~ s/$TEMPDIR/$TEMPURL/;
        $report_file =~ /([^\/]*$)/;
        $template->param(
            BODY => gen_body(
                report_file => $report_file,
                hsp_num     => $hsp_num,
                db_file     => $dbfile
            )
        );
    }
    else {
        $template->param( BODY => "Can't read or find file $report_file" );
    }
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    my $html;
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my %opts        = @_;
    my $report_file = $opts{report_file};
    my $db_file     = $opts{db_file};
    my $hsp_num     = $opts{hsp_num};
    my $hsps        = get_info_from_db(
        db_file     => $db_file,
        hsp_num     => $hsp_num,
        report_file => $report_file
    );
    my ( $qname, $sname ) = $hsps->[0]{hsp} =~ />?\(?(.*?)-(.*)\)?<?/;
    $qname =~ s/\(\(/\(/g;
    $sname =~ s/\(\(/\(/g;
    $qname =~ s/\)\)/\)/g;
    $sname =~ s/\)\)/\)/g;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'HSPView.tmpl' );
    $template->param( Query   => $qname );
    $template->param( Subject => $sname );
    my $qgap = $hsps->[0]{alignment} =~ tr/-/-/;
    my $sgap = $hsps->[1]{alignment} =~ tr/-/-/;

    $template->param( QGap => $qgap );
    $template->param( SGap => $sgap );
    my $qseq = $hsps->[0]{alignment};
    $qseq =~ s/-//g;
    my $sseq = $hsps->[1]{alignment};
    $sseq =~ s/-//g;
    $template->param( QLoc => $hsps->[0]{start} . "-" . $hsps->[0]{stop} );
    $template->param( SLoc => $hsps->[1]{start} . "-" . $hsps->[1]{stop} );
    my ($qmatch) = $hsps->[0]{match} =~ /(\d+)/;
    my ($smatch) = $hsps->[1]{match} =~ /(\d+)/;
    $template->param(
        QIdent => sprintf( "%.1f", $qmatch / length($qseq) * 100 ) . "%" )
      if $qseq;
    $template->param(
        SIdent => sprintf( "%.1f", $smatch / length($sseq) * 100 ) . "%" )
      if $sseq;
    $template->param(
        QSim => sprintf( "%.1f", $qmatch / length($qseq) * 100 ) . "%" )
      if $qseq;
    $template->param(
        SSim => sprintf( "%.1f", $smatch / length($sseq) * 100 ) . "%" )
      if $sseq;
    $template->param( QLen   => length($qseq) );
    $template->param( SLen   => length($sseq) );
    $template->param( QMis   => length($qseq) - $qmatch - $qgap );
    $template->param( SMis   => length($sseq) - $smatch - $sgap );
    $template->param( strand => $hsps->[0]{orientation} );
    $template->param( eval   => $hsps->[0]{eval} );
    $template->param( Match  => $hsps->[0]{match} );
    $template->param( Sim    => $hsps->[0]{match} );

    my $qseq_out = "<pre>" . seqwrap($qseq) . "</pre>";

    $template->param( qseq       => $qseq_out );
    $template->param( qseq_plain => $qseq );
    $template->param( sseq       => "<pre>" . seqwrap($sseq) . "</pre>" );
    $template->param( sseq_plain => $sseq );
    my @qln = split /\n/, seqwrap( $hsps->[0]{alignment}, 100 );
    my @sln = split /\n/, seqwrap( $hsps->[1]{alignment}, 100 );
    my @aln = split /\n/,
      seqwrap( get_alignment( $hsps->[0]{alignment}, $hsps->[1]{alignment} ),
        100 );
    my $align;

    for ( my $i = 0 ; $i < @aln ; $i++ ) {
        $align .= join( "\n", $qln[$i], $aln[$i], $sln[$i] ) . "\n\n";
    }
    $template->param( align => "<pre>" . $align . "</pre>" );
    if ( $hsps->[0]{dsid} =~ /^\d+$/ ) {
        my ( $start, $stop );
        $hsps->[0]{chr_start} =~ s/,//g;
        $hsps->[1]{chr_start} =~ s/,//g;
        $hsps->[0]{chr_stop}  =~ s/,//g;
        $hsps->[1]{chr_stop}  =~ s/,//g;
        $hsps->[0]{start}     =~ s/,//g;
        $hsps->[1]{start}     =~ s/,//g;
        $hsps->[0]{stop}      =~ s/,//g;
        $hsps->[1]{stop}      =~ s/,//g;
        $start = $hsps->[0]{start};
        $stop  = $hsps->[0]{stop};

        if ( $hsps->[0]{reverse_complement} ) {
            $start++;
            $stop++;
        }
        else {
            $start--;
            $stop--;
        }
        my $rc = $hsps->[0]{reverse_complement};
        my $seqview =
            "<a class = small href=SeqView.pl?dsid="
          . $hsps->[0]{dsid} . ";chr="
          . $hsps->[0]{chr}
          . ";start="
          . $start
          . ";stop="
          . $stop
          . ";rc=$rc>SeqView</a>";
        $template->param( qseqview => $seqview );
    }
    else {
        $template->param( qseqview => "<span class=small>N/A</small>" );
    }

    if ( $hsps->[1]{dsid} =~ /^\d+$/ ) {
        my ( $start, $stop );
        $start = $hsps->[1]{start};
        $stop  = $hsps->[1]{stop};

        if ( $hsps->[1]{reverse_complement} ) {
            $start++;
            $stop++;
        }
        else {
            $start--;
            $stop--;
        }
        my $rc = $hsps->[1]{reverse_complement};
        if ( $hsps->[1]{orientation} =~ /-/ ) {
            $rc = $rc ? 0 : 1;    #change orientation
        }
        my $seqview =
            "<a class = small href=SeqView.pl?dsid="
          . $hsps->[1]{dsid} . ";chr="
          . $hsps->[1]{chr}
          . ";start="
          . $start
          . ";stop="
          . $stop
          . ";rc=$rc>SeqView</a>";
        $template->param( sseqview => $seqview );
    }
    else {
        $template->param( sseqview => "<span class=small>N/A</small>" );
    }

    my $html;
    $html .= $template->output;
    return $html;
}

sub get_alignment {
    my @seq1 = split //, shift;
    my @seq2 = split //, shift;
    my $aln;
    for ( my $i = 0 ; $i < @seq1 ; $i++ ) {
        $aln .= $seq1[$i] eq $seq2[$i] ? "|" : " ";
    }
    return $aln;
}

sub seqwrap {
    my $seq     = shift;
    my $wrap    = shift || 50;
    my @chars   = split( //, $seq );
    my $newseq  = "";
    my $counter = 0;
    foreach my $char (@chars) {
        if ( $counter < $wrap ) {
            $newseq .= $char;
        }
        else {
            $newseq .= "\n" . $char;
            $counter = 0;
        }
        $counter++;
    }
    return ($newseq);
}

sub get_info_from_db {
    my %opts        = @_;
    my $db          = $opts{db_file};
    my $hsp_num     = $opts{hsp_num};
    my $report_file = $opts{report_file};
    my ($base) = $db =~ /^(.*?)\./;
    my ( $set1, $set2 ) = $report_file =~ /_(\d+)-(\d+)/;
    my $output;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$TEMPDIR/GEvo/$db", "", "" );
    my $statement =
qq{SELECT annotation, pair_id, id, image_id from image_data where name like "$hsp_num-%" and (image_id = $set1 or image_id = $set2)};
    my $sth = $dbh->prepare($statement);
    my $statement2 =
qq{SELECT dsid, chromosome, bpmin, bpmax, reverse_complement from image_info where id = ?};
    my $sth2 = $dbh->prepare($statement2);
    $sth->execute();
    my %hsps;
    my %ids;

    while ( my $item = $sth->fetchrow_arrayref ) {
        my %data;
        my $split;    # = "<.*?>";
        $split = "<br>"   if $item->[0] =~ /<br>/;
        $split = "\n"     if $item->[0] =~ /\n/;
        $split = "\\n"    if $item->[0] =~ /\\n/;
        $split = "<br\/>" if $item->[0] =~ /<br\/>/;
        $split = "&#10;"  if $item->[0] =~ /&#10;/;
        $split = "<tr>"   if $item->[0] =~ /<tr>/i;
        foreach my $set ( split /$split/i, $item->[0] ) {
            $set =~ s/<.*?>//g;
            next unless $set;
            my ( $tmp1, $tmp2 ) = split /:/, $set, 2;
            next unless $tmp2;
            $tmp2 =~ s/^\s+//;
            $tmp2 =~ s/\s+$//;
            $data{$tmp1} = $tmp2;
        }

        $data{HSP} =~ s/<.*?>//g;
        $data{HSP} =~ s/^\d+\s*//;
        my ( $start, $stop, $orientation ) =
          $data{Location} =~ /(.+?)-(.+?)\s+\(?(.*)\)?/;
        $orientation =~ s/\(|\)//g;
        $sth2->execute( $item->[3] );
        my ( $dsid, $chr, $bpmin, $bpmax, $rc ) = $sth2->fetchrow_array;
        $hsps{ $item->[2] } = {
            hsp                => $data{HSP},
            start              => $start,
            stop               => $stop,
            orientation        => $orientation,
            alignment          => $data{"Aligned Sequence"},
            match              => $data{Match},
            length             => $data{Length},
            identity           => $data{Identity},
            eval               => $data{E_val},
            location           => $data{Location},
            image_id           => $item->[3],
            dsid               => $dsid,
            chr_start          => $bpmin,
            chr_stop           => $bpmax,
            chr                => $chr,
            reverse_complement => $rc,
        };
        $ids{ $item->[1] }++;
        $ids{ $item->[2] }++;
    }
    my @hsps;

    #uberlame workaround ,but hey, it works
    foreach
      my $id ( sort { $hsps{$a}{image_id} <=> $hsps{$b}{image_id} } keys %hsps )
    {
        next unless $ids{$id} == 2;
        push @hsps, $hsps{$id};
    }
    return \@hsps;
}
