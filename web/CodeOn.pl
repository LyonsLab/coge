#! /usr/bin/perl -w

use strict;
use CoGeX;
use CoGeX::Result::Feature;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use CGI;
use CGI::Ajax;
use Data::Dumper;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use Benchmark;
use DBIxProfiler;
use CoGe::Accessory::genetic_code;
use Statistics::Basic::Mean;
use POSIX;
no warnings 'redefine';

use vars qw(
    $P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $PAGE_TITLE $TEMPDIR 
    $USER $DATE $BASEFILE $coge $cogeweb $FORM $COOKIE_NAME
);


$PAGE_TITLE = 'CodeOn';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$ENV{PATH} = $P->{COGEDIR};

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$TEMPDIR   = $P->{TEMPDIR};

$COOKIE_NAME = $P->{COOKIE_NAME};

# my ($cas_ticket) = $FORM->param('ticket');
# $USER = undef;
# ($USER) = CoGe::Accessory::Web->login_cas4(
#     cookie_name => $COOKIE_NAME,
#     ticket      => $cas_ticket,
#     coge        => $coge,
#     this_url    => $FORM->url()
# ) if ($cas_ticket);
# ($USER) = CoGe::Accessory::Web->get_user(
#     cookie_name => $COOKIE_NAME,
#     coge        => $coge
# ) unless $USER;

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $pj = new CGI::Ajax(
    gen_html => \&gen_html,
    get_orgs => \&get_orgs,
    go       => \&go,
);
$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );

#print $FORM->header,gen_html();

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'Coding Sequence Evolution',
                      PAGE_TITLE => 'CodeOn',
                      SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
                      HOME       => $P->{SERVER},
                      HELP       => 'CodeOn',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER       => $USER->display_name || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE       => $DATE );
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;

    return $html;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CodeOn.tmpl' );
    my $form = $FORM;
    $template->param( INITIALIZE => 1 );
    $template->param( ACCN => $form->param('accn') ) if $form->param('accn');
    $template->param( ANNO => $form->param('anno') ) if $form->param('anno');
    my @fids;
    @fids = map { split /::/, $_ } $form->param('fid');
    my @oids;
    @oids = map { split /::/, $_ } $form->param('oid');
    my @dsgids;
    @dsgids = map { split /::/, $_ } $form->param('dsgid');
    my $html;
    if ($USER->user_name eq "public") {
	    $html = "Due to the computationally intensive nature of CodeOn, you must be logged in to use this tool.";
    }	
    else {
        if ( @fids || @oids || @dsgids ) {
            my $res = go(
                fids   => \@fids,
                oids   => \@oids,
                dsgids => \@dsgids,
            );
            $template->param( RESULTS => $res );
        }
        $html .= $template->output;
    }
    return $html;
}

sub get_orgs {
    my %opts    = @_;
    my $type    = $opts{type};
    my $search  = $opts{search};
    my $id_only = $opts{id_only};
    my @db;
    my $count;
    if ( $type && $type eq "name" ) {
        @db =
          $coge->resultset("Organism")
          ->search( { name => { like => "%" . $search . "%" } } );
        $count = scalar @db;
    }
    elsif ( $type && $type eq "desc" ) {
        @db =
          $coge->resultset("Organism")
          ->search( { description => { like => "%" . $search . "%" } } );
        $count = scalar @db;
    }
    else {
        $count = $coge->resultset("Organism")->count();

        #       @db = $coge->resultset("Organism")->all;
    }
    return map { $_->id } @db if $id_only;

    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
        push @opts,
            "<OPTION value=\""
          . $item->id
          . "\" id=\"o"
          . $item->id . "\">"
          . $item->name
          . "</OPTION>";
    }
    my $html;
    $html .=
        qq{<FONT CLASS ="small" id="org_count">Organism count: }
      . $count
      . qq{</FONT>\n<BR>\n};
    if ( $search && !@opts ) {
        $html .= qq{<input type = hidden name="org_id" id="org_id"><br>};
        $html .= "No results";
        return $html;
    }
    unshift( @opts,
        "<OPTION value=\"all\" id=\"all\">All Listed Organisms</OPTION>" );
    my $size = scalar @opts;
    $size = 8 if $size > 8;
    $html .= qq{<SELECT id="org_id" SIZE="$size" MULTIPLE >\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
}

sub go {
    my %opts     = @_;
    my $accn     = $opts{accn};
    my $anno     = $opts{anno};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    my $oids     = $opts{oids};
    my $fids     = $opts{fids};
    my $dsgids   = $opts{dsgids};
    my $oid      = $opts{oid};        #scalar

    $fids   = [] unless defined $fids;
    $oids   = [] unless defined $oids;
    $dsgids = [] unless defined $dsgids;

    push @$oids, $oid if $oid;

    my ( $data, $feats, $dsgs ) = get_features(
        accn     => $accn,
        anno     => $anno,
        oids     => $oids,
        org_name => $org_name,
        org_desc => $org_desc,
        fids     => $fids,
        dsgids   => $dsgids
    );
    return $data unless ref($data) =~ /hash/i;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CodeOn.tmpl' );
    $template->param( RESULTS => 1 );
    my $aa_sort    = CoGe::Accessory::genetic_code->sort_aa_by_gc();
    my $table_head = "<th>"
      . join( "<th>",
        "GC% (feat count)",
        map    { $_ . "% (" . $data->{$_}{bin_count} . ")" }
          sort { $a <=> $b } keys %$data );
    $template->param( GC_HEAD => $table_head );
    my $max_aa = 0;
    my $min_aa = 100;

    foreach my $aa ( keys %$aa_sort ) {
        next if $aa eq "*";
        my (@tmp) =
          map  { $data->{$_}{data}{$aa} }
          sort { $data->{$b}{data}{$aa} <=> $data->{$a}{data}{$aa} }
          keys %$data;
        next unless defined $tmp[0];
        $max_aa = $tmp[0]  if $tmp[0] > $max_aa;
        $min_aa = $tmp[-1] if $tmp[-1] < $min_aa;
    }
    my @rows;
    foreach my $aa (
        sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
        keys %$aa_sort
      )
    {
        my @row;
        next if $aa eq "*";
        push @row,
          "<td>$aa (GC:" . sprintf( "%.0f", 100 * $aa_sort->{$aa} ) . "%)";
        foreach my $bin ( sort { $a <=> $b } keys %$data ) {
            my $aa_val = sprintf( "%.2f", 100 * $data->{$bin}{data}{$aa} );
            my $rel_val =
              ( $data->{$bin}{data}{$aa} - $min_aa ) / ( $max_aa - $min_aa );
            my $color = get_color( val => $rel_val );
            push @row,
                "<td style=\"color: black; background-color: rgb("
              . join( ",", @$color ) . ")\">"
              . $aa_val . "%";
        }
        push @rows, { RESULTS_ROW => join( "\n", @row ) };
    }
    my $featlistlink = "FeatList.pl?fid=" . join( "::", keys %$feats );
    my $dsg_info;
    if ( keys %$dsgs ) {
        $dsg_info = qq{<span class=small>};
        $dsg_info .= join(
            "<br>",
            map {
qq{<span class=link onclick="window.open('OrganismView.pl?dsgid=$_')">}
                  . $dsgs->{$_}->organism->name . " (v"
                  . $dsgs->{$_}->version
                  . qq{)</span>}
              } sort {
                $dsgs->{$a}->organism->name cmp $dsgs->{$b}->organism->name
              } keys %$dsgs
        );
        $dsg_info .= "</span>";
    }
    $template->param( FEAT_COUNT    => scalar( keys %$feats ) );
    $template->param( INFO          => \@rows );
    $template->param( FEATLISTLINK  => $featlistlink );
    $template->param( 'GENOME_INFO' => $dsg_info ) if $dsg_info;
    return $template->output;
}

sub get_features {
    my %opts     = @_;
    my $accn     = $opts{accn};
    my $anno     = $opts{anno};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    my $fids     = $opts{fids};       #array ref
    my $oids     = $opts{oids};       #array ref
    my $dsgids   = $opts{dsgids};     #array ref

    $fids   = [] unless defined $fids;
    $oids   = [] unless defined $oids;
    $dsgids = [] unless defined $dsgids;

    $org_name = undef if $org_name && $org_name =~ /search/i;
    $org_desc = undef if $org_desc && $org_desc =~ /search/i;
    my $weak_query = "Query needs to be better defined.";
    if (   !$accn
        && !$anno
        && !$fids
        && !$org_name
        && !$org_desc
        && !$oids
        && !$dsgids )
    {
        return $weak_query;
    }

    my $search = {};
    $search->{feature_type_id} = 3;
    $search->{'genome.organism_id'} = { IN => $oids } if $oids && @$oids;

    #    unless ($org_id)
    #      {
    $search->{'organism.name'} = { like => "%" . $org_name . "%" } if $org_name;
    $search->{'organism.description'} = { like => "%" . $org_desc . "%" }
      if $org_desc;

    #      }
    my $join = {
        join => [
            {
                feature => {
                    'dataset' =>
                      { 'dataset_connectors' => { 'genome' => 'organism' } }
                }
            }
        ]
    };

    #    push @{$join->{join}}, 'annotations' if $anno;#=>['annotation',]};
    #    push @{$join->{join}}, 'feature_names' if $accn;#=>['annotation',]};

    my %feats;
    my %orgs;
    my %dsgs;
    my $t1 = new Benchmark;
    if ($fids) {
        foreach my $fid (@$fids) {
            next unless $fid && $fid =~ /^\d+$/;
            $feats{$fid} = 1;

            #	    my $feat = $coge->resultset('Feature')->find($fid);
            #	    next unless $feat;
            #	    next unless $feat->feature_type_id == 3;
            #	    $feats{$feat->id} = $feat;
        }
    }
    if ($oids) {
        foreach my $oid (@$oids) {
            next unless $oid && $oid =~ /^\d+$/;
            $orgs{$oid} = 1;
        }
    }
    if ($dsgids) {
        foreach my $dsgid (@$dsgids) {
            next unless $dsgid && $dsgid =~ /^\d+$/;
            $dsgs{$dsgid} = 1;
        }
    }

    if ($accn) {
        map { $feats{ $_->feature->id } = $_->feature }
          $coge->resultset('FeatureName')->search( $search, $join )
          ->search_literal( 'MATCH(me.name) AGAINST (?)', $accn );
    }
    if ($anno) {
        map { $feats{ $_->feature->id } = $_->feature }
          $coge->resultset('FeatureAnnotation')->search( $search, $join )
          ->search_literal( 'MATCH(annotation) AGAINST (?)', $anno );
    }

    if ( $oids && @$oids && $oids->[0] eq "all" ) {
        my ( $otype, $search ) = ( "name", $org_name )
          if $org_name && $org_name ne "Search";
        ( $otype, $search ) = ( "desc", $org_desc )
          if $org_desc && $org_desc ne "Search";
        my @org_ids =
          get_orgs( id_only => 1, type => $otype, search => $search );
        map { $orgs{$_} = 1 } @org_ids;
    }

    foreach my $oid ( keys %orgs ) {
        my $org = $coge->resultset('Organism')->find($oid);
        next unless $org;
        foreach my $dsg ( $org->genomes ) {
            next unless $USER->has_access_to_genome($dsg);
            $dsgs{ $dsg->id } = $dsg;
        }
    }

    foreach my $dsgid ( keys %dsgs ) {
        my $dsg =
            $dsgs{$dsgid} eq "1"
          ? $coge->resultset('Genome')->find($dsgid)
          : $dsgs{$dsgid};
        next unless $dsg;
        $dsgs{$dsgid} = $dsg;
        foreach my $ds ( $dsg->datasets ) {
            map { $feats{ $_->id } = $_ }
              $ds->features( { feature_type_id => 3 } );
        }
    }

    foreach my $fid ( keys %feats ) {
        $feats{$fid} = $coge->resultset('Feature')->find($fid)
          if $feats{$fid} eq "1";
    }

    my %data;
    my $count = 0;
    my @aa;
    
    my $t2 = new Benchmark;

    my %codes;    #store genetic codes for reuse;
    foreach my $feat ( values %feats ) {
        my ($code) =
            $codes{ $feat->dataset->id }
          ? $codes{ $feat->dataset->id }
          : $feat->genetic_code;
        $codes{ $feat->dataset->id } = $code;
        my ( $gc, $at ) = $feat->gc_content();
        next unless $gc && $at;
        $gc *= 100;
        my $div = floor( $gc / 5 );
        my $mod = $gc % 5;

        #	print join ("\t", $gc, $div, $mod),"\n";
        $div++ if $mod > 2;
        my $aa_usage = $feat->aa_frequency( counts => 1, code => $code );
        @aa = keys %$aa_usage unless @aa;
        foreach my $aa ( keys %$aa_usage ) {
            push @{ $data{ 5 * $div }{$aa} }, $aa_usage->{$aa};
        }
        $count++;

        #	push @feats, {f=>$feat, gc=>$gc, at=>$at};
        #	last if $count > 50;
    }

    my $t3 = new Benchmark;

    my %return_data;
    foreach my $bin ( sort keys %data ) {
        $return_data{$bin}{bin_count} = scalar @{ $data{$bin}{"A"} };
        foreach my $aa (@aa) {
            my $ave = sprintf( "%.4f",
                Statistics::Basic::Mean->new( $data{$bin}{$aa} )->query );
            $return_data{$bin}{data}{$aa} = $ave;
        }
        $return_data{$bin}{data}{"*"} = 0 unless $return_data{$bin}{data}{"*"};
    }

    my $t4 = new Benchmark;

    foreach my $bin ( keys %return_data ) {
        my $hash  = $return_data{$bin}{data};
        my $total = 0;
        map { $total += $hash->{$_} } keys %$hash;
        unless ($total) {
            delete $return_data{$bin};
            next;
        }
        map { $hash->{$_} = $hash->{$_} / $total } keys %$hash;
        map { $total += $hash->{$_} } keys %$hash;
    }
    my $t5 = new Benchmark;

    my $query_time = timestr (timediff ($t2,$t1) );
    my $feat_time = timestr (timediff ($t3,$t2) );
    my $stats_time = timestr (timediff ($t4,$t3) );
    my $final_stats_time = timestr (timediff ($t5,$t4) );
    my $total_time = timestr( timediff ($t5, $t1) );
    print STDERR qq{

---CodeOn Benchmarks---
Time to do query:          $query_time
Time to process features:  $feat_time
Time for basic statistics: $stats_time
Time for final statistics: $final_stats_time
Total time:                $total_time

};
    return \%return_data, \%feats, \%dsgs;
}

sub get_color {
    my %opts = @_;
    my $val  = $opts{val};
    return [ 0, 0, 0 ] unless defined $val;
    return [ 125, 125, 125 ] if $val < 0;
    my @colors = (
        [ 255, 0,   0 ],      #red
        [ 255, 126, 0 ],      #orange
        [ 255, 255, 0 ],      #yellow
        [ 0,   255, 0 ],      # green
        [ 0,   255, 255 ],    # cyan
        [ 0,   0,   255 ],    # blue

        #                 [255,0,255], #magenta
        [ 126, 0, 126 ],      #purple
    );
    @colors = reverse @colors;
    my ( $index1, $index2 ) = (
        ( floor( ( scalar(@colors) - 1 ) * $val ) ),
        ceil( ( scalar(@colors) - 1 ) * $val )
    );

    my $color = [];
    my $step  = 1 / ( scalar(@colors) - 1 );
    my $scale = $index1 * $step;
    my $dist  = ( $val - $scale ) / ($step);
    for ( my $i = 0 ; $i <= 2 ; $i++ ) {
        my $diff = ( $colors[$index1][$i] - $colors[$index2][$i] ) * $dist;
        push @$color, sprintf( "%.0f", $colors[$index1][$i] - $diff );
    }
    return $color;
}
