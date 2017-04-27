#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;

use Benchmark;
use File::Path;
use DBI;
use POSIX;
use Digest::MD5 qw(md5_base64);
no warnings 'redefine';

use CoGe::Accessory::Web;
use CoGe::Accessory::genetic_code;
use CoGe::Accessory::Utils qw( commify );
use CoGeX;

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $DIAGSDIR $coge $COOKIE_NAME);
$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$DIAGSDIR  = $P->{DIAGSDIR};

# set this to 1 to print verbose messages to logs
$DEBUG = 0;

#$TEMPDIR = $P->{TEMPDIR}."SynSub";
#$TEMPURL = $P->{TEMPURL}."SynSub";
#mkpath ($TEMPDIR, 0,0777) unless -d $TEMPDIR;
$|    = 1;         # turn off buffering
$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM   = new CGI;
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

my $pj = new CGI::Ajax( gen_data => \&gen_data, );
$pj->JSDEBUG(0);
$pj->DEBUG(0);

#print $pj->build_html($FORM, \&gen_html);
print "Content-Type: text/html\n\n";
print gen_html($FORM);

sub gen_html {
    my ($body) = gen_body();

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => 'SynSub',
                      TITLE      => 'Synteny Substitution Matrix',
                      SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
                      HEAD       => qq{},
                      HOME       => $P->{SERVER},
                      HELP       => 'SynSub',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER       => $USER->display_name || '' );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE => $DATE );
    $template->param( BODY     => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    my $html = $template->output;
    return $html;
}

sub gen_body {
    my $form = shift || $FORM;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SynSub.tmpl' );

    my $dsgid1 = $form->param('dsgid1');
    my $dsgid2 = $form->param('dsgid2');
    my $db = $form->param('file');

    return "Need two dsgids to run." unless $dsgid1 && $dsgid2;
    my ($dsg1) = $coge->resultset('Genome')->find($dsgid1);
    my ($dsg2) = $coge->resultset('Genome')->find($dsgid2);

    my $color_type = $form->param('ct');    #color type for matrices output
    $color_type = "RYB" unless $color_type;

    #get gc content of genomes
    unless ( $dsg1 && $dsg2 ) {
        return
"<span class=alert>Problem generating dataset group objects for ids:  $dsgid1, $dsgid2.</span>";
    }
    my $org_name1 = $dsg1->organism->name;
    my $org_name2 = $dsg2->organism->name;
    ( $dsgid1, $dsg1, $org_name1, $dsgid2, $dsg2, $org_name2 ) =
      ( $dsgid2, $dsg2, $org_name2, $dsgid1, $dsg1, $org_name1 )
      if ( $org_name2 lt $org_name1 );
    my ( $dir1, $dir2 ) = sort ( $dsgid1, $dsgid2 );

    my $basedir  = "$DIAGSDIR/$dir1/$dir2";
    my $sqlite   = $basedir . "/" . $db;
    my (
        $freq_dna,               $org1_dna_counts,
        $org2_dna_counts,        $total_dna_counts,
        $freq_prot,              $org1_prot_counts,
        $org2_prot_counts,       $total_prot_counts,
        $org1_percent_gc_cds,    $org2_percent_gc_cds,
        $org1_percent_gc_wobble, $org2_percent_gc_wobble
    ) = get_counts( sqlite => $sqlite );
    my ( $dna_matrix1, $dna_matrix2 ) = gen_matrix(
        freq        => $freq_dna,
        org1_counts => $org1_dna_counts,
        org2_counts => $org2_dna_counts
    );
    my ( $prot_matrix1, $prot_matrix2 ) = gen_matrix(
        freq        => $freq_prot,
        org1_counts => $org1_prot_counts,
        org2_counts => $org2_prot_counts
    );

    my ( $org1_gc, $org1_length ) = get_gc_dsg($dsg1);
    my ( $org2_gc, $org2_length ) = get_gc_dsg($dsg2);
    my $temp_name1 = $org_name1;
    my $temp_name2 = $org_name2;
    $temp_name1 =~ s/(^\S+\s+\S+)\s/$1<br>/;
    $temp_name2 =~ s/(^\S+\s+\S+)\s/$1<br>/;
    my $org_info;
    $template->param( 'DSGID1'         => $dsgid1 );
    $template->param( 'DSGID2'         => $dsgid2 );
    $template->param( 'ORG1_NAME'      => $temp_name1 );
    $template->param( 'ORG2_NAME'      => $temp_name2 );
    $template->param( 'ORG1_LENGTH'    => $org1_length );
    $template->param( 'ORG2_LENGTH'    => $org2_length );
    $template->param( 'ORG1_GC'        => $org1_gc );
    $template->param( 'ORG2_GC'        => $org2_gc );
    $template->param( 'ORG1_GC_CDS'    => $org1_percent_gc_cds );
    $template->param( 'ORG2_GC_CDS'    => $org2_percent_gc_cds );
    $template->param( 'ORG1_GC_WOBBLE' => $org1_percent_gc_wobble );
    $template->param( 'ORG2_GC_WOBBLE' => $org2_percent_gc_wobble );
    my $synmap_link = "SynMap.pl?dsgid1=$dsgid1;dsgid2=$dsgid2;ks=1;autogo=1";
    $synmap_link =
qq{<a href='$synmap_link' class='coge-button' style='color: #000000' target=_new_synmap>SynMap</a>};
    my $synsub_link = "SynSub.pl?dsgid1=$dsgid1;dsgid2=$dsgid2;ct=";
    $synsub_link =
        qq{<a href='$synsub_link} . 'rain'
      . qq{' class='coge-button' style='color: #000000' target=_new_synmap>Rainbown SynSub</a>}
      . qq{<a href='$synsub_link} . 'RYB'
      . qq{' class='coge-button' style='color: #000000' target=_new_synmap>RYB SynSub</a>};

    $template->param( LINKS => $synmap_link . $synsub_link );
    my $html;
    $html .=
      "<table class='ui-widget ui-corner-all ui-widget-content small'><tr>";
    $html .= "<td>Scores normalized to protein frequences of $org_name1<br>";
    $html .= gen_matrix_output_html(
        matrix      => $prot_matrix1,
        log         => 0,
        type        => "protein",
        org1        => $org_name1,
        org2        => $org_name2,
        org1_counts => $org1_prot_counts,
        color_type  => $color_type
    );
    $html .= "<td width=1px>";
    $html .= "<td>Scores normalized to protein frequences of $org_name2<br>";
    $html .= gen_matrix_output_html(
        matrix      => $prot_matrix2,
        log         => 0,
        type        => "protein",
        org1        => $org_name1,
        org2        => $org_name2,
        org1_counts => $org2_prot_counts,
        color_type  => $color_type
    );
    $html .= "</table>";
    $html .= "<hr>";
    $html .= "Scores normalized to codon frequences of $org_name1<br>";
    $html .= gen_matrix_output_html(
        matrix      => $dna_matrix1,
        log         => 0,
        type        => "codon",
        org1        => $org_name1,
        org2        => $org_name2,
        org1_counts => $org1_dna_counts,
        color_type  => $color_type
    );
    $html .= "<hr>";
    $html .= "Scores normalized to codon frequences of $org_name2<br>";
    $html .= gen_matrix_output_html(
        matrix      => $dna_matrix2,
        log         => 0,
        type        => "codon",
        org1        => $org_name1,
        org2        => $org_name2,
        org1_counts => $org2_dna_counts,
        color_type  => $color_type
    );
    $html .= "<hr>";
    $template->param( MATRIX => $html );
    return $template->output;
}

sub get_counts {
    my %opts   = @_;
    my $sqlite = $opts{sqlite};

    my %freq_dna;
    my $org1_dna_counts  = {};
    my $org2_dna_counts  = {};
    my $total_dna_counts = {};

    my %freq_prot;
    my $org1_prot_counts  = {};
    my $org2_prot_counts  = {};
    my $total_prot_counts = {};

    my $org1_cds_gc       = 0;
    my $org2_cds_gc       = 0;
    my $org1_wobble_gc    = 0;
    my $org2_wobble_gc    = 0;
    my $org1_cds_total    = 0;
    my $org2_cds_total    = 0;
    my $org1_wobble_total = 0;
    my $org2_wobble_total = 0;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$sqlite", "", "" );
    my $select =
qq{SELECT protein_align_1, protein_align_2, DNA_align_1, DNA_align_2 FROM ks_data};
    my $sth = $dbh->prepare($select);
    $sth->execute;
    while ( my $item = $sth->fetchrow_arrayref() ) {
        next unless $item->[2];

        my $dna_length = length( $item->[2] );
        my ( $length, $gc );
        ( $length, $gc ) = get_gc_seq( seq => $item->[2] );
        $org1_cds_total += $length;
        $org1_cds_gc    += $gc;
        ( $length, $gc ) = get_gc_seq( seq => $item->[3] );
        $org2_cds_total += $length;
        $org2_cds_gc    += $gc;

        #wobbles
        ( $length, $gc ) = get_wobble_gc_seq( seq => $item->[2] );
        $org1_wobble_total += $length;
        $org1_wobble_gc    += $gc;
        ( $length, $gc ) = get_wobble_gc_seq( seq => $item->[3] );
        $org2_wobble_total += $length;
        $org2_wobble_gc    += $gc;

        for ( my $i = 0 ; $i < $dna_length ; $i += 3 )    #getting codons
        {
            my $chr1 = uc( substr( $item->[2], $i, 3 ) );
            my $chr2 = uc( substr( $item->[3], $i, 3 ) );
            next if $chr1 =~ /-/;
            next if $chr2 =~ /-/;
            next if $chr1 =~ /[^ATCG]/;
            next if $chr2 =~ /[^ATCG]/;
            $freq_dna{$chr1}
              {$chr2}++;    #increment frequency count of align characters
        }

        my $prot_length = length( $item->[0] );
        for ( my $i = 0 ; $i < $prot_length ; $i += 1 )    #getting amino acids
        {
            my $chr1 = uc( substr( $item->[0], $i, 1 ) );
            my $chr2 = uc( substr( $item->[1], $i, 1 ) );
            next if $chr1 =~ /-/;
            next if $chr2 =~ /-/;
            $freq_prot{$chr1}
              {$chr2}++;    #increment frequency count of align characters
        }
    }

    foreach my $c1 ( keys %freq_dna ) {
        if (
            length($c1) !=
            3 )             #sometimes, something goes terribly, terribly wrong
        {
            delete $freq_dna{$c1};
            next;
        }
        foreach my $c2 ( keys %{ $freq_dna{$c1} } ) {
            if (
                length($c2) !=
                3 )         #sometimes, something goes terribly, terribly wrong
            {
                delete $freq_dna{$c1}{$c2};
                next;
            }

            $org1_dna_counts->{$c1}  += $freq_dna{$c1}{$c2};
            $org2_dna_counts->{$c2}  += $freq_dna{$c1}{$c2};
            $total_dna_counts->{$c1} += $freq_dna{$c1}{$c2};
            $total_dna_counts->{$c2} += $freq_dna{$c1}{$c2};
        }
    }

    foreach my $c1 ( keys %freq_prot ) {
        foreach my $c2 ( keys %{ $freq_prot{$c1} } ) {
            $org1_prot_counts->{$c1}  += $freq_prot{$c1}{$c2};
            $org2_prot_counts->{$c2}  += $freq_prot{$c1}{$c2};
            $total_prot_counts->{$c1} += $freq_prot{$c1}{$c2};
            $total_prot_counts->{$c2} += $freq_prot{$c1}{$c2};
        }
    }
    my $org1_percent_gc_cds =
      sprintf( "%.2f", 100 * $org1_cds_gc / $org1_cds_total );
    my $org2_percent_gc_cds =
      sprintf( "%.2f", 100 * $org2_cds_gc / $org2_cds_total );
    my $org1_percent_gc_wobble =
      sprintf( "%.2f", 100 * $org1_wobble_gc / $org1_wobble_total );
    my $org2_percent_gc_wobble =
      sprintf( "%.2f", 100 * $org2_wobble_gc / $org2_wobble_total );

    return \%freq_dna, $org1_dna_counts, $org2_dna_counts, $total_dna_counts,
      \%freq_prot, $org1_prot_counts, $org2_prot_counts, $total_prot_counts,
      $org1_percent_gc_cds, $org2_percent_gc_cds, $org1_percent_gc_wobble,
      $org2_percent_gc_wobble;
}

sub gen_matrix {
    my %opts = @_;

    my $freq        = $opts{freq};
    my $org1_counts = $opts{org1_counts};
    my $org2_counts = $opts{org2_counts};

#log-odds score based on Henikoff and Henikoff 1992 paper:  Amino acid substitution matrices from protein blocks
#calculate the total number of aligned characters
    my $org1_total = 0;
    my $org2_total = 0;
    map { $org1_total += $_ } values %$org1_counts;
    map { $org2_total += $_ } values %$org2_counts;

    #org1_total and org2_total should be the same value!
    my $total_pairs = $org1_total + $org2_total;
    my %q
      ; #count of $chr1 changing to $chr2 divided by the total number of aligned characters $total_pairs
    my @chrs = keys %$freq;
    foreach my $c1 (@chrs) {
        foreach my $c2 (@chrs) {
            $q{$c1}{$c2} = $freq->{$c1}{$c2} / $org1_total if $freq->{$c1}{$c2};
        }
    }

    #calcuate the total fequency of $chr1
    my %p1;
    my %p2;
    foreach my $c1 (@chrs) {
        foreach my $c2 (@chrs) {
            $p1{$c1} += $q{$c1}{$c2} if $q{$c1}{$c2};
            $p2{$c2} += $q{$c1}{$c2} if $q{$c1}{$c2};
        }
    }

    my %matrix1;
    my %matrix2;

    foreach my $c1 (@chrs) {
        foreach my $c2 (@chrs) {

            #$q is the frequency of $c1 changing to $c2
            my $q = $q{$c1}{$c2} if $q{$c1}{$c2};
            next unless $q;
            my $denom =
              ( $p1{$c1} * $p1{$c2} );    #combined frequency of $c1 * $c2
            my $val = sprintf( "%.2f", 2 * log( $q / $denom ) ) if $denom;
            $val = 0 if $val eq "-0";
            $matrix1{$c1}{$c2} = $val;
            $denom = ( $p2{$c1} * $p2{$c2} )
              if defined $p2{$c1}
                  && defined $p2{$c2};    #combined frequency of $c1 * $c2
            $val = sprintf( "%.2f", 2 * log( $q / $denom ) ) if $denom;
            $val = 0 if $val eq "-0";
            $matrix2{$c1}{$c2} = $val;
        }
    }

    #    print Dumper \%matrix1, \%matrix2;
    #    return \%q1, \%q2;
    return \%matrix1, \%matrix2;          #, \%q, \%p;
}

sub gen_matrix_output_html {
    my %opts        = @_;
    my $data        = $opts{matrix};
    my $type        = $opts{type};       #c for codons, p for protein, d for dna
    my $log         = $opts{log};        #flag to log transform data values
    my $org1        = $opts{org1};
    my $org2        = $opts{org2};
    my $org1_counts = $opts{org1_counts};
    my $org2_counts = $opts{org2_counts};
    my $color_type  = $opts{color_type};
    my $code = CoGe::Accessory::genetic_code->code;
    $code = $code->{code};
    my $html;
    my @order;

    if ( $type =~ /c/i ) {
        my @dna = qw(A T C G);
        my @codons;
        foreach my $c1 (@dna) {
            foreach my $c2 (@dna) {
                foreach my $c3 (@dna) {
                    push @codons, $c1 . $c2 . $c3;
                }
            }
        }
        @order = sort {
            sort_nt1( substr( $a, 1, 1 ) ) <=> sort_nt1( substr( $b, 1, 1 ) )
              || sort_nt2( substr( $a, 0, 1 ) )
              <=> sort_nt2( substr( $b, 0, 1 ) )
              || sort_nt3( substr( $a, 2, 1 ) )
              <=> sort_nt3( substr( $b, 2, 1 ) )
        } @codons;
    }
    else {
        my $aa_sort = CoGe::Accessory::genetic_code->sort_aa_by_gc();
        @order =
          sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
          keys %$aa_sort;
    }

#    my %check = map {$_=>1} @order;  #used to remove data that will not be displayed.

#some entries will not be used and need to be skipped.  Examples are amino acid "X" and nucleotide "N" for unknowns
    my ( $min, $max );
    foreach my $c1 (@order) {
        foreach my $c2 (@order) {
            my $val = $data->{$c1}{$c2};
            next unless defined $val;
            $val = $log ? log10($val) : $val;
            $max = $val unless defined $max;
            $max = $val if $val > $max;
            $min = $val unless defined $min;
            $min = $val if $val < $min;
        }
    }
    my $range = $max - $min;

    #fill in items without values
    foreach my $c1 (@order) {
        foreach my $c2 (@order) {
            $data->{$c1}{$c2} = "*" unless defined $data->{$c1}{$c2};
        }
    }

    $html .=
qq{<table class='ui-widget ui-corner-all ui-widget-content small' style='border-collapse: collapse;'>};
    $html .= "<tr><th>";
    my $col_count = 0;
    foreach my $item (@order) {
        $col_count++;
        $html .= "<th>$item";
        $html .= "<br>" . "(" . $code->{$item} . ")" if $code->{$item};
        $html .= "<td bgcolor=grey width=1px>"
          unless $col_count % 16 || $type =~ /p/;
    }
    $html .= "<th>Total";
    $html .= "<tr>";
    my $row_count = 0;
    foreach my $aa1 (@order) {
        $col_count = 0;
        $row_count++;
        $html .= "<th>$aa1";
        $html .= "(" . $code->{$aa1} . ")" if $code->{$aa1};
        my %vals;
        map { $vals{$_}++ } map { $data->{$aa1}{$_} } keys %{ $data->{$aa1} };
        my $total = 0;
        foreach my $aa2 (@order)    #	    foreach my $aa2 (sort keys %$data)
        {
            $col_count++;
            my $val = $data->{$aa1}{$aa2};
            next unless $range;
            my $relative_val = $val =~ /\d/ ? ( $val - $min ) / $range : $val;
            $val = sprintf( "%.1f", $val ) if $val =~ /\d/;
            my $color = get_color( val => $relative_val, type => $color_type );
            my $color_str = join( ",", @$color );
            my $font_color = "#000000";
            $font_color = "#AAAAAA"
              if $color->[0] <= 85 && $color->[1] <= 85 && $color->[2] <= 150;
            $font_color = "#AAAAAA" if $color->[0] <= 100 && $color->[1] <= 50;
            $html .=
"<td style=\"background-color: rgb($color_str); color: $font_color\">"
              . $val;    #." ".$code->{$aa1}."-".$code->{$aa2};
            $html .= "<td bgcolor=grey width=1px>"
              unless $col_count % 16 || $type =~ /p/;
        }
        $html .= "<td align=right>";
        $html .= $org1_counts->{$aa1} ? commify( $org1_counts->{$aa1} ) : 0;
        $html .= "<tr bgcolor=grey width=1px><td colspan=80>"
          unless $row_count % 16 || $type =~ /p/;
        $html .= "<tr>";
    }

    #    $html .= "<th>Total";
    #    $col_count =0;
    #    foreach my $item (@order)
    #      {
    #	$col_count ++;
    #	$html .= "<td>";
    #	$html .= $org2_counts->{$item} ? $org2_counts->{$item} : 0;;
    #	$html .= "<td bgcolor=grey width=5px>" unless $col_count %16;
    #      }
    $html .= "</table>";
    $html .= "x-axis: " . $org2 . "<br>";
    $html .= "y-axis: " . $org1 . "<br>";
    return $html;
}

sub get_color {
    my %opts = @_;
    my $val  = $opts{val};
    my $type = $opts{type};
    return [ 0, 0, 0 ] unless defined $val;
    return [ 200, 200, 200 ] if $val !~ /\d/;
    my @rainbow = (
        [ 255, 0,   0 ],      #red
        [ 255, 255, 0 ],      #yellow
        [ 0,   255, 0 ],      # green
        [ 0,   255, 255 ],    # cyan
        [ 220, 0,   220 ],    #magenta
        [ 0,   0,   255 ],    # blue
    );
    my @red_yellow_blue = (
        [ 255, 0,   0 ],      #red
        [ 220, 220, 20 ],     #yellow
        [ 0,   0,   150 ],    # blue
    );
    my @colors;
    @colors = $type =~ /rain/ ? @rainbow : @red_yellow_blue;
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

sub sort_nt1 {
    my $chr = uc(shift);

    $chr = substr( $chr, -1, 1 ) if length($chr) > 1;
    my $val = 0;
    if ( $chr eq "C" ) {
        $val = 1;
    }
    elsif ( $chr eq "T" ) {
        $val = 2;
    }
    elsif ( $chr eq "A" ) {
        $val = 3;
    }
    return $val;
}

sub sort_nt2 {
    my $chr = uc(shift);

    $chr = substr( $chr, -1, 1 ) if length($chr) > 1;
    my $val = 0;
    if ( $chr eq "G" ) {
        $val = 1;
    }
    elsif ( $chr eq "A" ) {
        $val = 2;
    }
    elsif ( $chr eq "T" ) {
        $val = 3;
    }
    return $val;
}

sub sort_nt3 {
    my $chr = uc(shift);

    $chr = substr( $chr, -1, 1 ) if length($chr) > 1;
    my $val = 0;
    if ( $chr eq "G" ) {
        $val = 1;
    }
    elsif ( $chr eq "T" ) {
        $val = 2;
    }
    elsif ( $chr eq "C" ) {
        $val = 3;
    }
    return $val;
}

sub get_gc_dsg {
    my $dsg      = shift;
    my $length   = 0;
    my $gc_total = 0;
    foreach my $chr ( $dsg->chromosomes ) {
        my ( $gc, $at, $n, $x ) = $dsg->percent_gc( count => 1, chr => $chr );
        $gc_total += $gc;
        $length   += ( $gc + $at + $n + $x );
    }
    my $percent_gc = sprintf( "%.2f", 100 * $gc_total / $length );
    return $percent_gc, commify($length);
}

sub get_gc_seq {
    my %opts = @_;
    my $seq  = $opts{seq};
    $seq = uc($seq);
    $seq =~ s/-//g;    #remove gaps;
    my ($gc) = $seq =~ tr/GC/GC/;
    return ( length($seq), $gc );
}

sub get_wobble_gc_seq {
    my %opts = @_;
    my $seq  = $opts{seq};
    $seq = uc($seq);
    $seq =~ s/-//g;    #remove gaps;
    my $length = length($seq);
    my $wobbles;
    for ( my $i = 0 ; $i < $length ; $i += 3 ) {
        next if $i + 2 > $length;
        my $chr = substr( $seq, $i + 2, 1 );
        $wobbles .= $chr;
    }
    return ( get_gc_seq( seq => $wobbles ) );
}

1;
