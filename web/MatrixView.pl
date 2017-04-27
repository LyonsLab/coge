#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGeX;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;

#use CoGeX;
use Benchmark;
use CoGe::Accessory::Web;
use CoGe::Accessory::genetic_code;
use Statistics::Basic::Mean;
use Digest::MD5 qw(md5_base64);
use POSIX;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $TEMPDIR $TEMPURL $MATRIXDIR $USER $FORM $coge $connstr $COOKIE_NAME);

# set this to 1 to print verbose messages to logs
$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$DEBUG     = 0;
$TEMPDIR   = $P->{TEMPDIR};
$TEMPURL   = $P->{TEMPURL};
$MATRIXDIR = $P->{BLASTMATRIX} . "aa/";
$|         = 1;                                      # turn off buffering
$DATE      = sprintf(
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
($USER) = CoGe::Accessory::Web->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

my $pj = new CGI::Ajax( gen_data => \&gen_data, );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html( $FORM, \&gen_html );

#print "Content-Type: text/html\n\n";print gen_html($FORM);

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'Sequence Alignment Matrix View',
                      PAGE_TITLE => 'MatrixView',
                      HOME       => $P->{SERVER},
                      HELP       => 'MatrixView',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      HEAD       => qq{},
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
    my $form = shift || $FORM;
    my $matrix = $form->param('matrix');
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'MatrixView.tmpl' );
    my $matrices = process_dir( matrix => $matrix );
    $template->param( MATRICES => $matrices );
    my $data = gen_data( file => $matrix );
    $template->param( DATA => $data );
    return $template->output;

}

sub gen_data {
    my %opts        = @_;
    my $form        = $opts{form} || $FORM;
    my $code_layout = $opts{code_layout} || 0;
    my $bin_size    = $opts{bin_size} || 5;
    my $file        = $opts{file};
    return unless $file;
    my $code = CoGe::Accessory::genetic_code->code;
    $code = $code->{code};
    my $html = "<table>";
    my ($data) = process_file( file => $file );
    my ( $max, $min );

    foreach my $c1 ( keys %$data ) {
        foreach my $val ( values %{ $data->{$c1} } ) {
            next if $val == -50;
            next if $val > 25;
            $max = $val unless defined $max;
            $max = $val if $val > $max;
            $min = $val unless defined $min;
            $min = $val if $val < $min;
        }
    }
    my $range = $max - $min;
    if ( $data->{P} ) {
        my $aa_sort = CoGe::Accessory::genetic_code->sort_aa_by_gc();
        $html .= "<tr><th><th>"
          . join( "<th>",
            sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
              keys %$aa_sort );
        $html .= "<th>Total:";
        $html .= "<tr>";
        foreach my $aa1 (
            sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
            keys %$aa_sort
          )
        {
            $html .= "<th>$aa1";
            my $total = 0;
            foreach my $aa2 (
                sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
                keys %$aa_sort
              )
            {
                my $val = $data->{$aa1}{$aa2};
                $total += $val;
                my $relative_val = ( $val - $min ) / $range;
                my $color        = get_color( val => $relative_val );
                my $color_str    = join( ",", @$color );
                $html .=
                  "<td style=\"background-color: rgb($color_str)\">" . $val;

            }
            $html .= "<td>$total<tr>";
        }
    }
    else {
        my @order = sort {
            sort_nt1( substr( $a, 1, 1 ) ) <=> sort_nt1( substr( $b, 1, 1 ) )
              || sort_nt2( substr( $a, 0, 1 ) )
              <=> sort_nt2( substr( $b, 0, 1 ) )
              || sort_nt3( substr( $a, 2, 1 ) )
              <=> sort_nt3( substr( $b, 2, 1 ) )
        } keys %$data;
        $html .= "<tr><th>";
        my $col_count = 0;
        foreach my $item (@order) {
            $col_count++;
            $html .= "<th>$item";
            $html .= "<td bgcolor=grey width=5px>" unless $col_count % 16;
        }
        $html .= "<th>Total";
        $html .= "<tr>";
        my $row_count = 0;
        foreach my $aa1 (@order) {
            $col_count = 0;
            $row_count++;
            $html .= "<th>$aa1";
            my $total = 0;
            foreach my $aa2 (@order) {
                $col_count++;
                my $val = $data->{$aa1}{$aa2};
                $total += $val;
                my $relative_val = ( $val - $min ) / $range;
                my $color        = get_color( val => $relative_val );
                my $color_str    = join( ",", @$color );
                $html .=
                  "<td style=\"background-color: rgb($color_str)\">" . $val;
                $html .= "<td bgcolor=grey width=5px>" unless $col_count % 16;
            }
            $html .= "<td>$total";
            $html .= "<tr bgcolor=grey width=5px><td colspan=80>"
              unless $row_count % 16;
            $html .= "<tr>";
        }
    }
    $html .= "</table>";
    return $html;
}

sub process_dir {
    my %opts   = @_;
    my $matrix = $opts{matrix};
    my %matrices;
    opendir( DIR, $MATRIXDIR ) || warn $!;
    while ( my $file = readdir(DIR) ) {
        next if $file =~ /^\./;
        $matrices{"$file"} = uc($file)
          unless $matrices{ uc($file) } || $matrices{ lc($file) };
    }
    closedir DIR;
    my $html = "<select id=matrix>";
    $html .= join( "\n",
        map  { "<option value=\"" . $_ . "\">" . $matrices{$_} . "</option>" }
        sort { $matrices{$a} cmp $matrices{$b} } keys %matrices )
      . "\n";
    $html =~ s/(value=\"$matrix\")/selected $1/i if $matrix;
    $html .= "</select>";
    return $html;
}

sub process_file {
    my %opts = @_;
    my $bin_size = $opts{bin_size} || 5;
    $bin_size = 1   if $bin_size < 1;
    $bin_size = 100 if $bin_size > 100;
    my $skip = $opts{skip} || [];    #[qw(mitochond chloroplast virus phage)];
    my $keep = $opts{keep} || [];    #[qw(mitochondr)];
    my $file = $opts{file};
    $file = $MATRIXDIR . "/" . $file unless $file =~ /$MATRIXDIR/i;
    my %data;
    my @head;
    my $max = 0;
    my $min = 10000;
    open( IN, $file );

    while (<IN>) {
        chomp;
        next if /^#/;
        my @line = split /\s+/;
        unless ( $line[0] ) {
            shift @line;
            @head = @line;
            next;
        }
        my $aa = shift @line;
        for ( my $i = 0 ; $i <= $#line ; $i++ ) {
            $data{$aa}{ $head[$i] } = $line[$i];
            $max = $line[$i] if $line[$i] > $max;
            next if $aa eq "*";
            next if $head[$i] eq "*";
            $min = $line[$i] if $line[$i] < $min;
        }
    }
    close IN;
    return \%data, $max, $min;
}

sub get_color {
    my %opts = @_;
    my $val  = $opts{val};
    return [ 0, 0, 0 ] unless defined $val;
    return [ 125, 125, 125 ] if $val < 0;
    return [ 125, 125, 125 ] if $val > 1;
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

1;
