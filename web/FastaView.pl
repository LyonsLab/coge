#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGeX;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use File::Basename;
use Text::Wrap qw($columns &wrap);
use JSON::XS;

no warnings 'redefine';

use vars qw($P $TEMPDIR $TEMPURL $FORM $USER $LINK $coge $PAGE_TITLE $PAGE_NAME);

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$PAGE_TITLE = 'FastaView';
$PAGE_NAME  = "$PAGE_TITLE.pl";
$ENV{PATH}  = $P->{COGEDIR};
$TEMPDIR    = $P->{TEMPDIR} . "$PAGE_TITLE/";
$TEMPURL    = $P->{TEMPURL} . "$PAGE_TITLE/";

my %ajax = CoGe::Accessory::Web::ajax_func();

my %FUNCTIONS = (
    gen_html => \&gen_html,
    get_seqs => \&get_seqs,
    gen_file => \&gen_file,
    test => \&test,
    %ajax,
);

#CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );
my $pj = new CGI::Ajax(%FUNCTIONS);
$pj->js_encode_function('escape');

if ( $FORM->param('jquery_ajax') ) {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};

    #print STDERR Dumper \%args;
    if ( $fname and defined $FUNCTIONS{$fname} ) {
        if ( $args{args} ) {
            my @args_list = split( /,/, $args{args} );
            print $FORM->header, $FUNCTIONS{$fname}->(@args_list);
        }
        else {
            print $FORM->header, $FUNCTIONS{$fname}->(%args);
        }
    }
}
elsif ( $FORM->param('text') ) {
    my $header =
      "Content-disposition: attachement; filename=CoGe_";    #test.gff\n\n";
    $header .= int( rand(100000) );
    $header .= ".faa\n\n";

    print $header;
    print gen_html();
}
else {
    print $pj->build_html( $FORM, \&gen_html );
    #print $FORM->header,gen_html();
}

sub gen_html {
    my $html;
    my $form       = $FORM;
    my $prot       = $form->param('prot');
    my $text       = $form->param('text');
    my $name_only  = $form->param('no');
    my $id_only    = $form->param('io');
    my $upstream   = $form->param('up') || 0;
    my $downstream = $form->param('down') || 0;
    my $textbox    = $text ? 0 : 1;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );

    #       $template->param(TITLE=>'Fasta Viewer');
    $template->param( PAGE_TITLE => 'FastaView',
				      TITLE      => 'FastaView',
    				  PAGE_LINK  => $LINK,
    				  HOME       => $P->{SERVER},
                      HELP       => 'FastaView',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER => $USER->display_name || '' );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    my @fids;
    push @fids, $form->param('featid') if $form->param('featid');
    push @fids, $form->param('fid')    if $form->param('fid');
    my $gstid = $form->param('gstid') if $form->param('gstid');

    my ( $seqs, $seq_count, $feat_count, $warning ) = get_seqs(
        prot       => $prot,
        fids       => \@fids,
        textbox    => $textbox,
        gstid      => $gstid,
        name_only  => $name_only,
        id_only    => $id_only,
        upstream   => $upstream,
        downstream => $downstream
    );

    my $json = get_json(
        prot       => $prot,
        fids       => \@fids,
        textbox    => $textbox,
        gstid      => $gstid,
        name_only  => $name_only,
        id_only    => $id_only,
        upstream   => $upstream,
        downstream => $downstream
    );

    if ($text) {
        return $seqs;
    }
    $template->param(
        BODY => gen_body(
            fids       => \@fids,
            seqs       => $seqs,
            json       => $json,
            seq_count  => $seq_count,
            feat_count => $feat_count,
            gstid      => $gstid,
            prot       => $prot,
            up         => $upstream,
            down       => $downstream,
            message    => $warning,
        )
    );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );

    $html .= $template->output;
    return $html;
}

sub get_seqs {
    my %opts       = @_;
    my $fids       = $opts{fids};
    my $prot       = $opts{prot};
    my $textbox    = $opts{textbox};
    my $name_only  = $opts{name_only};
    my $id_only    = $opts{id_only};
    my $gstid      = $opts{gstid};
    my $upstream   = $opts{upstream};
    my $downstream = $opts{downstream};
    my @fids       = ref($fids) =~ /array/i ? @$fids : split( /,/, $fids );
    my %seen       = ();
    @fids = grep { !$seen{$_}++ } @fids;

    my $seqs;
    my $seq_count = 0;
    my $fid_count = 0;
    foreach my $item (@fids) {
        foreach my $featid ( split( /,/, $item ) ) {
            $fid_count++;
            my ( $fid, $gstidt );
            if ( $featid =~ /_/ ) {
                ( $fid, $gstidt ) = split /_/, $featid;
            }
            else {
                ( $fid, $gstidt ) = ( $featid, $gstid );
            }
            my ($feat) = $coge->resultset('Feature')->find($fid);
            unless ($feat) {
                $seqs .= ">Not found: $featid\n";
                next;
            }
            my ($dsg) = $feat->dataset->genomes;
            if ( !$USER->has_access_to_genome($dsg) ) {
                $seqs .= ">Restricted: $featid\n";
                next;
            }
            my $tmp = $feat->fasta(
                col        => 100,
                prot       => $prot,
                name_only  => $name_only,
                fid_only   => $id_only,
                gstid      => $gstidt,
                upstream   => $upstream,
                downstream => $downstream
            );
            $seq_count += $tmp =~ tr/>/>/;
            $seqs .= $tmp;
        }
    }
    my $warning;
    %seen = ();
    while ( $seqs =~ /(>.*\n)/g ) {
        $warning = "Warning: Duplicate sequence names" if $seen{$1};
        $seen{$1} = 1;
    }
    $seqs =
qq{<textarea id=seq_text name=seq_text class="ui-widget-content ui-corner-all backbox" ondblclick="this.select();" style="height: 400px; width: 800px; overflow: auto;">$seqs</textarea>}
      if $textbox;
    return $seqs, $seq_count, $fid_count, $warning;
}

sub get_json { # mdb added 2/28/14 for genfam integration
    my %opts       = @_;
    my $fids       = $opts{fids};
    my $prot       = $opts{prot};
    my $gstid      = $opts{gstid};
    my $upstream   = $opts{upstream};
    my $downstream = $opts{downstream};
    my @fids       = ref($fids) =~ /array/i ? @$fids : split( /,/, $fids );
    my %seen       = ();
    @fids = grep { !$seen{$_}++ } @fids;

    my @obj;
    foreach my $item (@fids) {
        foreach my $featid ( split( /,/, $item ) ) {
            my ( $fid, $gstidt );
            if ( $featid =~ /_/ ) {
                ( $fid, $gstidt ) = split /_/, $featid;
            }
            else {
                ( $fid, $gstidt ) = ( $featid, $gstid );
            }
            my ($feat) = $coge->resultset('Feature')->find($fid);
            next unless ($feat);

            my ($dsg) = $feat->dataset->genomes;
            next if ( !$USER->has_access_to_genome($dsg) );

            my $tmp = $feat->fasta_object(
                col        => 100,
                prot       => $prot,
                gstid      => $gstidt,
                upstream   => $upstream,
                downstream => $downstream
            );
            push @obj, $tmp;
        }
    }

    return encode_json({fasta => \@obj});
}

sub gen_file {
    my %opts       = @_;
    my $fids       = $opts{fids};
    my $prot       = $opts{prot};
    my $textbox    = $opts{textbox};
    my $name_only  = $opts{name_only};
    my $id_only    = $opts{id_only};
    my $gstid      = $opts{gstid};
    my $upstream   = $opts{upstream};
    my $downstream = $opts{downstream};
    my @fids       = ref($fids) =~ /array/i ? @$fids : split(/,/, $fids);

    my ($seqs)     = get_seqs(
        prot       => $prot,
        fids       => \@fids,
        gstid      => $gstid,
        name_only  => $name_only,
        id_only    => $id_only,
        upstream   => $upstream,
        downstream => $downstream
    );

    my $file = $TEMPDIR . "Seqs_" . int( rand(100000000) ) . ".faa";
    eval {
        open( OUT, ">" . $file );
        print OUT $seqs;
        close OUT;
    };

    return $P->{SERVER} . "/tmp/$PAGE_TITLE/" . basename($file);
}

sub gen_body {
    my %opts       = @_;
    my $seqs       = $opts{seqs};
    my $json       = $opts{json};
    my $seq_count  = $opts{seq_count};
    my $feat_count = $opts{feat_count};
    my $message    = $opts{message};
    my $fids       = $opts{fids};
    my $gstid      = $opts{gstid} || 1;
    my $prot       = $opts{prot} || 0;
    my $up         = $opts{up} || 0;
    my $down       = $opts{down} || 0;
    $fids = join( ",", @$fids ) if ref($fids) =~ /array/i;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'FastaView.tmpl' );
    $template->param( BOTTOM_BUTTONS => 1 );
    $template->param( SEQ            => $seqs ) if $seqs;
    $template->param( FASTA_JSON     => $json ) if $json;
    $template->param( SEQ_COUNT      => $seq_count ) if defined $seq_count;
    $template->param( FEAT_COUNT     => $feat_count ) if defined $feat_count;
    $template->param( WARNING        => $message ) if defined $message;
    $template->param( FIDS =>
qq{<input type=hidden id=fids value=$fids><input type=hidden id=gstid value=$gstid>}
    );
    $template->param( PROT   => $prot );
    $template->param( UP     => $up );
    $template->param( DOWN   => $down );
    $template->param( genfam => $P->{GENFAMURL});
    return $template->output;
}
