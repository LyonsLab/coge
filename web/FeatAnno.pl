#! /usr/bin/perl -w
use strict;

use CGI;
use CGI::Carp 'fatalsToBrowser';
use Data::Dumper;
use DBIxProfiler;
use Digest::MD5 qw(md5_base64);

use CoGeX;
use CoGe::Accessory::Web;
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG
  $TEMPDIR $TEMPURL $USER $FORM $FID $DS $DSG $CHR $LOC $ORG $TYPE
  $VERSION $START $STOP $NAME_ONLY $coge $GSTID $COOKIE_NAME);
$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$TEMPDIR   = $P->{TEMPDIR};
$TEMPURL   = $P->{TEMPURL};

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$|     = 1;    # turn off buffering

$FORM = new CGI;
$FID  = $FORM->param('fid');
$DS   = $FORM->param('ds');     # || 61;
$DSG  = $FORM->param('dsg');    # || 61;
$CHR  = $FORM->param('chr');    # || 7;
$LOC =
     $FORM->param('loc')
  || $FORM->param('pos')
  || $FORM->param('x');         # || 6049802;
$LOC       = 0 unless $LOC;
$START     = $FORM->param('start');
$START     = $LOC unless defined $START;
$STOP      = $FORM->param('stop');
$STOP      = $START unless defined $STOP;
$ORG       = $FORM->param('org') || $FORM->param('organism');
$VERSION   = $FORM->param('version') || $FORM->param('ver');
$NAME_ONLY = $FORM->param('name_only') || 0;
$GSTID  = $FORM->param('gstid');    #genomic_sequence_type_id
$TYPE   = $FORM->param('type');     # mdb added 6/18/13
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

print "Content-Type: text/html\n\n",
      '<link rel="stylesheet" type="text/css" href="css/coge.css" />'; # mdb added 5/9/13 for JBrowse
my $rhtml = gen_html(
    featid    => $FID,
    start     => $START,
    stop      => $STOP,
    chr       => $CHR,
    ds        => $DS,
    dsg       => $DSG,
    org       => $ORG,
    version   => $VERSION,
    name_only => $NAME_ONLY,
    gstid     => $GSTID
) if $START > 0 || $FID;

if ( $START && $STOP ) {
    if ( $START == $STOP ) {
        $rhtml =
"<font class=title3>Position:</font> <font class=data>$START</font><br><hr>"
          . $rhtml;
    }
    else {
        $rhtml =
"<font class=title3>Position:</font> <font class=data>$START-$STOP</font><br><hr>"
          . $rhtml;
    }
}
$rhtml = "No annotations" unless $rhtml;
print $rhtml;

sub gen_html {
    my %args   = @_;
    my $featid = $args{featid};
    my $start  = $args{start};
    my $stop   = $args{stop};
    $stop = $start unless $stop;
    my $chr       = $args{chr};
    my $ds        = $args{ds};
    my $dsg       = $args{dsg};
    my $version   = $args{version};
    my $name_only = $args{name_only};
    my $gstid     = $args{gstid};
    ($ds)  = $coge->resultset('Dataset')->resolve($ds) if $ds;
    ($dsg) = $coge->resultset('Genome')->resolve($dsg) if $dsg;
    my @ds;
    push @ds, $ds            if $ds;
    push @ds, $dsg->datasets if $dsg;

    if ( $ds
      ) #there can be additional information about a chromosome for a particular version of the organism that is not in the same data_information.  Let's go find the organism_id and version for the specified data information item.
    {
        $version = $ds->version unless $version;
        ($chr) = $ds->get_chromosomes() unless $chr;
    }
    if ($dsg) {
        $version = $dsg->version unless $version;
        ($chr) = $dsg->chromosomes() unless ($chr);
    }
    my @all_feats;
    foreach my $item (@ds) {
        push @all_feats,
          $coge->get_features_in_region(
            dataset_id => $item->id,
            chr        => $chr,
            start      => $start,
            stop       => $stop,
          ) if ( $chr && $start && $stop );
    }
    push @all_feats, $coge->resultset('Feature')->find($featid) if $featid;
    return " " unless @all_feats;

    my @feats;
    foreach my $feat (@all_feats) {
        my $type_name = $feat->type->name;
        next if $type_name =~ /^source$/;
        next if $type_name =~ /^chromosome$/;
        next if ( defined $TYPE and $type_name ne $TYPE );   # mdb added 6/18/13
        push @feats, $feat;
    }

    my $html;
    if ($name_only) {
        $html = "<table>";
        my $class = "even";
        foreach my $feat ( sort { $a->type->name cmp $b->type->name } @feats ) {
            $html .=
                "<tr class='small $class'><td class='title5'>"
              . $feat->type->name
              . ": </td><td >"
              . (
                join(
                    ", ",
                    map {
                        "<span onclick=\"window.open('FeatView.pl?accn=$_;fid="
                          . $feat->id
                          . "');\" class='link'>$_</span>"
                      } $feat->names
                )
              ) . "</td></tr>";
            $class = $class eq "even" ? "odd" : "even";
        }
        $html .= "</table>";
        return $html;
    }
    foreach my $feat ( sort { $a->type->name cmp $b->type->name } @feats ) {
        my $color;
        if ( $feat->type->name eq "CDS" ) {
            $color = "#DDFFDD";
        }
        elsif ( $feat->type->name eq "mRNA" ) {
            $color = "#DDDDFF";
        }
        elsif ( $feat->type->name =~ "rna" ) {
            $color = "#DDDDDD";
        }
        elsif ( $feat->type->name =~ "gene" ) {
            $color = "#FFDDDD";
        }
        else {
            $color = "#FFDDBB";
        }
        $html .= '<div style="background-color:';
        $html .= $color;
        $html .= ';width:100%"><h4>';
        $html .= $feat->type->name;
        $html .= '</h4>';
        $html .= $feat->annotation_pretty_print_html( gstid => $gstid, P => $P);
        $html .= '</div><hr>';
    }
    return $html;
}
