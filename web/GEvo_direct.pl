#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CGI;
use CGI::Ajax;
use Image::Size;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use DBI;
use File::Path;
use CoGe::Accessory::Web;
no warnings 'redefine';

delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw ($P $FORM $USER $TEMPDIR $TEMPURL $LINK $PAGE_TITLE $PAGE_NAME $coge);
$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};

$FORM    = new CGI;
$TEMPDIR = $P->{TEMPDIR} . "GEvo";
$TEMPURL = $P->{TEMPURL} . "GEvo";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;

$PAGE_TITLE = "GEvo_direct";
$PAGE_NAME  = "$PAGE_TITLE.pl";

( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my $pj = new CGI::Ajax(
    gen_html => \&gen_html,
    update   => \&update,
);

$pj->js_encode_function('escape');

#print $pj->build_html($FORM, \&gen_html);

print $FORM->header, gen_html();

sub gen_html {
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( BODY => gen_body() );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER       => $name,
                      HELP       => "/wiki/index.php?title=GEvo_direct",
                      TITLE      => "GEvo direct:  reviewing past results.", 
                      PAGE_TITLE => "GEvo direct" );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    
    return $template->output;
}

sub gen_body {
    my $form         = $FORM;
    my $name         = $form->param('name');
    my $cmd          = "/usr/bin/svnversion " . $P->{COGEDIR} . "gobe/flash";
    my $gobe_version = `$cmd`;
    $gobe_version =~ s/\n//g;

    my %files;
    my $tiny;

    open( CMD, "/bin/ls $TEMPDIR/$name" . "* |" );
    while (<CMD>) {
        my $touch = "/usr/bin/touch $_";
        my $x;
        ( $x, $touch ) = CoGe::Accessory::Web::check_taint($touch);
        `$touch`;
        foreach ( split /\n/ ) {
            if (/\.anno/) {
                push @{ $files{anno} }, $_;
            }
            elsif (/\.faa/) {
                push @{ $files{faa} }, $_;
            }
            elsif (/\.png/) {
                push @{ $files{png} }, $_;
            }
            elsif (/\.log/) {
                open( IN, $_ );
                while ( my $line = <IN> ) {
                    if ( $line =~ /tiny url: (.*)/i ) {
                        $tiny = $1;
                        last;
                    }
                }
                close IN;

                push @{ $files{log} }, $_;
            }
            elsif (/\.sqlite/) {
                push @{ $files{sqlite} }, $_;
            }
            else {
                push @{ $files{report} }, $_;
            }
        }
    }
    close CMD;

    #print STDERR Dumper \%files;
    my $h       = 0;
    my $w       = 0;
    my $seq_num = 0;
    foreach my $img ( @{ $files{png} } ) {
        my ( $x, $y ) = imgsize($img);
        $h += $y;
        $w = $x;
        $seq_num++;
    }
    $w += 5;
    my $html;
    my $dbname = $files{sqlite}[0];
    $html .= get_db_stuff($dbname);
    unless ( $html =~ /does not exist/ )    #check to see if results exists
    {
        $html .= qq{<table>};
        $html .= qq{<tr valign=top><td class = small>Alignment reports};
        my $i = 1;
        foreach my $report ( @{ $files{report} } ) {
            $report =~ s/$TEMPDIR/$TEMPURL/;
            $html .=
"<div><font class=small><A HREF=\"$report\" target=_new>View alignment output $i</A></font></DIV>\n";
            $i++;
        }

        $html .= qq{<td class = small>Fasta files};
        $i = 1;
        foreach my $item ( @{ $files{faa} } ) {
            $item =~ s/$TEMPDIR/$TEMPURL/;
            $html .=
"<div><font class=small><A HREF=\"$item\" target=_new>Fasta file $i</A></font></DIV>\n";
            $i++;
        }

        $html .=
qq{<td class = small><a href = "http://baboon.math.berkeley.edu/mavid/gaf.html">GAF</a> annotation files};
        $i = 1;
        foreach my $item ( @{ $files{anno} } ) {
            $item =~ s/$TEMPDIR/$TEMPURL/;
            $html .=
"<div><font class=small><A HREF=\"$item\" target=_new>Annotation file $i</A></font></DIV>\n";
            $i++;
        }
        $html .= qq{<td class = small>SQLite db};
        $dbname =~ s/$TEMPDIR/$TEMPURL/;
        $html .=
"<div class=small><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
        $html .= qq{<td class = small>Log File};
        my $logfile = $files{log}[0];
        $logfile =~ s/$TEMPDIR/$TEMPURL/;
        $html .=
          "<div class=small><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
        $html .= qq{<td class = small>GEvo Link<div class=small>};
        $html .= qq{<a href=$tiny target=_new>$tiny<br>} if $tiny;
        $html .= qq{(See log file for full link)</a></div>};
        $html .= qq{</table>};
    }
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GEvo_direct.tmpl' );
    $template->param( 'STUFF'        => $html );
    $template->param( 'WIDTH'        => $w );
    $template->param( 'HEIGHT'       => $h );
    $template->param( 'GOBE_VERSION' => $gobe_version );
    $template->param( 'SEQ_NUM'      => $seq_num );
    $template->param( 'BASE_NAME'    => $name );
    return $template->output;
}

sub get_db_stuff {
    my $dbname = shift;
    my $html;
    unless ( -r $dbname ) {
        $html .=
"Database file $dbname does not exist.  CoGe only stores the results from GEvo's analysis for ~24 hours, and your prior results may have been deleted.  For long term \"storage\" of your results, please copy the GEvo link, as a tinyurl, provided below the results or view the log file (a link to which is also provided below the results) which has the full GEvo url.";
        return $html;
    }

    $html .= qq{
<Div class="topItem">Change image display order:</DIV>
<div class="dropMenu">
<table>
 <tr>
  <th>Order<th>Image Name
};
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbname", "", "" );
    my $query = qq{
select * from image_info order by display_id asc
};
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $count = 1;

    #    my @stuff;
    my $display = '$(\'#display';
    my $image   = '$(\'#image';

    while ( my $data = $sth->fetchrow_arrayref ) {
        my $id         = $data->[0];
        my $display_id = $data->[1];
        my $name       = $data->[3];
        $html .=
qq{<tr><td><input type="text" id="display-$count" value="$count" size=2><input type=hidden id="image-$count" value=$id><td>$name};

  #	push @stuff, qq{'args__'+$display-$count').val()+'_'+$image-$count').val()};

        $count++;
    }
    $count--;
    $html .= qq{</table>};
    $html .= qq{<input type=hidden id="num_seqs" value="$count">};

    #    my $params = join (", ",@stuff);
    #    print STDERR $params;
    $dbname =~ s/$TEMPDIR\/?//;
    $html .= qq{<input type=hidden id="sqlite" value=$dbname>};

    $html .= qq{
<input type=button value="Update Image" onclick=prepare_update()>

};
    $html .= "</div>";
    $sth->finish();
    $sth = undef;
    $dbh->disconnect();
    return $html;
}

sub update {
    my $dbname = shift @_;
    my %data;
    my $count = 1;
    foreach my $item (
        sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
        map { [ split /-/ ] } split /:/,
        shift @_
      )
    {
        $data{$count} = $item->[1];
        $count++;
    }
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$TEMPDIR/$dbname", "", "" );
    while ( my ( $k, $v ) = each %data ) {
        my $query = qq{
UPDATE image_info set display_id = $k where id = $v;
};
        print STDERR $query unless $dbh->do($query);
    }
    $dbh->disconnect();
}
