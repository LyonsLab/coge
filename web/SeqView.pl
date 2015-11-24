#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use POSIX;
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );

use vars qw($P $PAGE_NAME $PAGE_TITLE $FORM $USER $coge $LINK);

$PAGE_TITLE = 'SeqView';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my $pj = new CGI::Ajax(
    gen_html           => \&gen_html,
    get_seq            => \&get_seq,
    gen_title          => \&gen_title,
    find_feats         => \&find_feats,
    parse_url          => \&parse_url,
    generate_feat_info => \&generate_feat_info,
    generate_gc_info   => \&generate_gc_info,
);
$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );

#print $FORM->header, gen_html();

sub gen_html {
    my $html;
    unless ($USER) {
        $html = login();
    }
    else {
        my $form = $FORM;
        my $rc   = $form->param('rc');
        my $pro;
        my ($title) = gen_title( protein => $pro, rc => $rc );
        my $template =
          HTML::Template->new(
            filename => $P->{TMPLDIR} . 'generic_page.tmpl' );

        $template->param( PAGE_TITLE => 'SeqView',
        				  PAGE_LINK  => $LINK,
        				  HOME       => $P->{SERVER},
                          HELP       => 'SeqView',
                          WIKI_URL   => $P->{WIKI_URL} || '',
                          USER => $USER->display_name || '' );
        $template->param( BODY       => gen_body() );
        $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
        $template->param( ADMIN_ONLY => $USER->is_admin );
        $template->param( CAS_URL    => $P->{CAS_URL} || '' );
        $html .= $template->output;
    }
    return $html;
}

sub gen_body {
    my $form   = $FORM;
    my $featid = $form->param('featid') || $form->param('fid') || 0;
    my $gstid  = $form->param('gstid') if $form->param('gstid');
    ( $featid, $gstid ) = split( /_/, $featid ) if ( $featid =~ /_/ );

    my $chr        = $form->param('chr');
    my $dsid       = $form->param('dsid');
    my $dsgid      = $form->param('dsgid');
    my $feat_name  = $form->param('featname');
    my $rc         = $form->param('rc');
    my $pro        = $form->param('pro');
    my $upstream   = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $locations  = $form->param('locations');
    my $rel        = $form->param('rel')
      || 0
      ; #relative position of the feature -- don't adjust based on which strand the feature is located.  important when linking to seqview from places such as GEvo's Get Sequence
    my $start = $form->param('start');
    $start =~ s/,//g  if $start;
    $start =~ s/\.//g if $start;
    my $stop = $form->param('stop');
    $stop =~ s/,//g  if $stop;
    $stop =~ s/\.//g if $stop;
    $stop = $start unless $stop;
    ( $start, $stop ) = ( $stop, $start ) if $start && $stop && $start > $stop;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SeqView.tmpl' );
    $template->param( RC       => $rc );
    $template->param( REL      => $rel );
    $template->param( JS       => 1 );
    $template->param( SEQ_BOX  => 1 );
    $template->param( ADDITION => 1 );
    $template->param( GSTID    => $gstid );
    $template->param( DSID     => $dsid );
    $template->param( DSGID    => $dsgid );
    $template->param( CHR      => $chr );

    if ($featid) {
        my ($feat) = $coge->resultset('Feature')->find($featid);
        $dsid = $feat->dataset_id;
        $chr  = $feat->chromosome;

        $template->param(
            FEAT_START => $feat->start,
            FEAT_STOP  => $feat->stop,
            FEATID     => $featid,
            FEATNAME   => $feat_name,
            PROTEIN    => 'Protein Sequence',
            SIXFRAME   => 0,
            UPSTREAM   => "Add 5': ",
            UPVALUE    => $upstream,
            DOWNSTREAM => "Add 3': ",
            DOWNVALUE  => $downstream,
            FEATURE    => 1
        );
        $start = $feat->start;
        $stop  = $feat->stop;
    }
    elsif ($locations) {
        $template->param(
            PROTEIN    => 'Six Frame Translation',
            SIXFRAME   => 0,
            UPSTREAM   => "Add 5': ",
            UPVALUE    => $upstream,
            DOWNSTREAM => "Add 3': ",
            DOWNVALUE  => $downstream,
            LOCATIONS  => $locations
        );
    }
    else {
        $template->param(
            PROTEIN    => 'Six Frame Translation',
            SIXFRAME   => 1,
            UPSTREAM   => "Start: ",
            UPVALUE    => $start,
            DOWNSTREAM => "Stop: ",
            DOWNVALUE  => $stop,
            ADD_EXTRA  => 1,
            ADDUP      => $upstream,
            ADDDOWN    => $downstream
        );

    }

    if ($rc) {
        $start -= $downstream;
        $stop += $upstream;
    }
    else {
        $start -= $upstream;
        $stop += $downstream;
    }

    my ( $link, $types ) = find_feats(
        dsid  => $dsid,
        start => $start,
        stop  => $stop,
        chr   => $chr,
        gstid => $gstid,
        dsgid => $dsgid
    );

    $template->param( FEATLISTLINK   => $link );
    $template->param( FEAT_TYPE_LIST => $types );

    return $template->output;
}

sub check_strand {
    my %opts   = @_;
    my $strand = $opts{'strand'} || 1;
    my $rc     = $opts{'rc'} || 0;
    if ( $rc == 1 ) {
        if ( $strand =~ /-/ ) {
            $strand = "1";
        }
        else {
            $strand = "-1";
        }
    }
    elsif ( $strand =~ /-/ ) {
        $strand =~ s/^\-$/-1/;
    }
    else {
        $strand =~ s/^\+$/1/;
    }
    return $strand;
}

sub get_seq {
    my %opts       = @_;
    my $add_to_seq = $opts{'add'};
    my $featid     = $opts{'featid'} || 0;
    $featid = 0 if $featid eq "undefined";    #javascript funkiness
    my $pro = $opts{'pro'};

    #my $pro = 1;
    my $rc         = $opts{'rc'} || 0;
    my $chr        = $opts{'chr'};
    my $dsid       = $opts{'dsid'};
    my $dsgid      = $opts{'dsgid'};
    my $feat_name  = $opts{'featname'};
    my $upstream   = $opts{'upstream'};
    my $downstream = $opts{'downstream'};
    my $locations  = $opts{'locations'};
    my $start      = $opts{'start'};
    my $stop       = $opts{'stop'};
    my $wrap       = $opts{'wrap'} || 0;
    my $gstid      = $opts{gstid};
    my $rel        = $opts{rel} || 0;
    $wrap = 0 if $wrap =~ /undefined/;

    if ($add_to_seq) {
        $start = $upstream   if $upstream;
        $stop  = $downstream if $downstream;
    }
    else {
        $start -= $upstream;
        $stop += $downstream;
    }
    my $strand;
    my $seq;
    my $fasta;
    my $col = $wrap ? 80 : 0;

    if ($featid) {
        my $feat = $coge->resultset('Feature')->find($featid);
        my ($dsg) = $feat->dataset->genomes;
        return "Restricted Access" unless $USER->has_access_to_genome($dsg);

#	return "Restricted Access" if $feat->dataset->restricted && !$USER->has_access_to_dataset($feat->dataset);
        ( $fasta, $seq ) =
          ref($feat) =~ /Feature/i
          ? $feat->fasta(
            prot       => $pro,
            rc         => $rc,
            upstream   => $upstream,
            downstream => $downstream,
            col        => $col,
            sep        => 1,
            gstid      => $gstid,
            rel        => $rel,
          )
          : ">Unable to retrieve Feature object for id: $featid\n";

#	$seq = $rc ? color(seq=>$seq, upstream=>$downstream, downstream=>$upstream) : color(seq=>$seq, upstream=>$upstream, downstream=>$downstream);
        $fasta = $fasta . "\n" . $seq . "\n";
    }
    elsif ($dsid) {
        my $ds = $coge->resultset('Dataset')->find($dsid);
        my ($dsg) = $ds->genomes;
        return "Restricted Access" unless $USER->has_access_to_genome($dsg);
        $fasta =
          ref($ds) =~ /dataset/i
          ? $ds->fasta(
            start => $start,
            stop  => $stop,
            chr   => $chr,
            prot  => $pro,
            rc    => $rc,
            col   => $col,
            gstid => $gstid,
          )
          : ">Unable to retrieve dataset object for id: $dsid";
    }
    elsif ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        return "Unable to find genome for $dsgid" unless $dsg;
        return "Restricted Access" unless $USER->has_access_to_genome($dsg);
        $fasta =
          ref($dsg) =~ /genome/i
          ? $dsg->fasta(
            start => $start,
            stop  => $stop,
            chr   => $chr,
            prot  => $pro,
            rc    => $rc,
            col   => $col,
          )
          : ">Unable to retrieve dataset group object for id: $dsgid";
    }
    elsif ($locations) {
        my %genomes;
        foreach ( split( ',', $locations ) ) {
            my ( $gid, $chr, $start, $stop ) = split( ':', $_ );
            next if ( $gid   =~ /\D/ );
            next if ( $start =~ /\D/ );
            next if ( $stop  =~ /\D/ );
            ( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
            $start -= $upstream;
            $stop += $downstream;
            if ( not defined $genomes{$gid} ) {
                $genomes{$gid} = $coge->resultset('Genome')->find($gid);
            }
            my $genome = $genomes{$gid};
            return "Unable to find genome for $gid" unless $genome;
            return "Restricted Access" unless $USER->has_access_to_genome($genome);
            $fasta .=
              ref($genome) =~ /genome/i
              ? $genome->fasta(
                start => $start,
                stop  => $stop,
                chr   => $chr,
                prot  => $pro,
                rc    => $rc,
                col   => $col,
              )
              : ">Unable to retrieve dataset group object for id: $gid";
        }
    }
    else {
        $fasta = qq{
>Unable to create sequence.  Options:
};
        $fasta .= Dumper \%opts;
    }
    return $fasta;
}

sub color {
    my %opts = @_;
    my $seq  = $opts{'seq'};

    #       my $rc = $opts{'rc'};
    my $upstream   = $opts{'upstream'};
    my $downstream = $opts{'downstream'};
    $upstream   = 0 if $upstream < 0;
    $downstream = 0 if $downstream < 0;
    my $up;
    my $down;
    my $main;
    my $nl1;
    $nl1 = 0;
    $up = substr( $seq, 0, $upstream );
    while ( $up =~ /\n/g ) { $nl1++; }
    my $check = substr( $seq, $upstream, $nl1 );

    $nl1++ if $check =~ /\n/;
    $upstream += $nl1;
    $up = substr( $seq, 0, $upstream );
    my $nl2 = 0;
    $down = substr( $seq, ( ( length $seq ) - ($downstream) ), length $seq );
    while ( $down =~ /\n/g ) { $nl2++; }
    $check = substr( $seq, ( ( length $seq ) - ( $downstream + $nl2 ) ), $nl2 );

    $nl2++ if $check =~ /\n/;
    $downstream += $nl2;
    $down = substr( $seq, ( ( length $seq ) - ($downstream) ), $downstream );
    $up   = lc($up);
    $down = lc($down);
    $main =
      substr( $seq, $upstream,
        ( ( ( length $seq ) ) - ( $downstream + $upstream ) ) );
    $main = uc($main);
    $seq = join( "", $up, $main, $down );
    return $seq;
}

sub gen_title {
    my %opts     = @_;
    my $rc       = $opts{'rc'} || 0;
    my $pro      = $opts{'pro'};
    my $sixframe = $opts{sixframe};
    my $title;
    if ($pro) {
        $title = $sixframe ? "Six Frame Translation" : "Protein Sequence";
    }
    else {
        $title = $rc ? "Reverse Complement" : "DNA Sequence";
    }
    return $title;
}

sub find_feats {
    my %opts  = @_;
    my $start = $opts{'start'};
    my $stop  = $opts{'stop'};
    my $chr   = $opts{'chr'};
    my $dsid  = $opts{'dsid'};
    my $gstid = $opts{'gstid'};
    my $dsgid = $opts{dsgid};
    my @dsids;
    push @dsids, $dsid if $dsid;

    return unless ($chr);

    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        return unless $dsg;
        push @dsids, map { $_->id } $dsg->datasets;
        $gstid = $dsg->type->id;
    }
    my $link =
qq{<span class='ui-button ui-corner-all' " onClick="featlist('FeatList.pl?};
    my %type;
    $link .=
        "start=$start;stop=$stop;chr=$chr;dsid=$dsid;dsgid=$dsgid;gstid=$gstid"
      . qq{')">Extract Features:};
    foreach my $ft (
        $coge->resultset('FeatureType')->search(
            {
                "features.dataset_id" => [@dsids],
                "features.chromosome" => $chr
            },
            {
                join   => "features",
                select => [ { "distinct" => "me.feature_type_id" }, "name" ],
                as => [ "feature_type_id", "name" ],
            }
        )
      )
    {
        $type{ $ft->name } = $ft->id;
    }
    $type{All} = 0;

    my $type = qq{<SELECT ID="feature_type">};
    $type .= join( "\n",
        map { "<OPTION value=" . $type{$_} . ">" . $_ . "</option>" }
        sort keys %type )
      . "\n";
    $type .= "</select>";
    return $link, $type;
}

sub generate_feat_info {
    my $featid = shift;
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless ( ref($feat) =~ /Feature/i ) {
        return "Unable to retrieve Feature object for id: $featid";
    }

#    my $html = qq{<a href="#" onClick="\$('#feature_info').slideToggle(pageObj.speed);" style="float: right;"><img src='/CoGe/picts/delete.png' width='16' height='16' border='0'></a>};
    my $html = $feat->annotation_pretty_print_html();
    return $html;
}

sub generate_gc_info {
    my $seq      = shift;
    my $seq_type = shift;
    return "Cannot Calculate GC content of Protein Sequence" if $seq_type;
    $seq =~ s/>.*?\n//;
    $seq =~ s/\n//g;
    my $length = length($seq);
    return "No sequence" unless $length;
    my $gc = $seq =~ tr/GCgc/GCgc/;
    my $at = $seq =~ tr/ATat/ATat/;
    my $n  = $seq =~ tr/Nn/Nn/;
    my $x  = $seq =~ tr/Xx/Xx/;
    my $pgc = sprintf( "%.2f", $gc / $length * 100 );
    my $pat = sprintf( "%.2f", $at / $length * 100 );
    my $pn  = sprintf( "%.2f", $n / $length * 100 );
    my $px  = sprintf( "%.2f", $x / $length * 100 );
    my $total_content = qq{<table class=small>\n};
    $total_content .= "<tr align=right><th>NT<th>Count<th>Percent";
    $total_content .=
        "<tr align=right><td>GC:<td align=right>"
      . commify($gc) . "<td>("
      . $pgc . "%)";
    $total_content .=
        "<tr align=right><td>AT:<td align=right>"
      . commify($at) . "<td>("
      . $pat . "%)";
    $total_content .=
        "<tr align=right><td>N:<td align=right>"
      . commify($n) . "<td>("
      . $pn . "%)";
    $total_content .=
        "<tr align=right><td>X:<td align=right>"
      . commify($x) . "<td>("
      . $px . "%)";
    $total_content .=
      "<tr align=right><td>Total:<td align=right>" . commify($length);
    $total_content .= qq{</table>\n};
    return $total_content;
}
