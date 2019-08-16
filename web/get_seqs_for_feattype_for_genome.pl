#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use Data::Dumper;
use Text::Wrap;
use CGI;
use IO::Compress::Gzip;
use File::Path;

$|++;

use vars qw($FORM $P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $coge $FASTADIR $URL $DIR $USER $COOKIE_NAME $DB $USER $config $LINK);

$FORM = new CGI;
$P    = CoGe::Accessory::Web::get_defaults();

# EHL added 8.16.19
( $DB, $USER, $config, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => "get_seqs_for_feattype_for_genome"
);
#end add

$FASTADIR = $P->{FASTADIR};
$DIR      = $P->{COGEDIR};
$URL      = $P->{URL};
mkpath( $FASTADIR, 1, 0777 );
unless (-r $FASTADIR) {
	print $FORM->header;
	print "Access denied to $FASTADIR\n";
	exit;
}

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
unless ($coge) {
	print $FORM->header;
	print "Failed to connect to database\n";
	exit;
}

#$COOKIE_NAME = $P->{COOKIE_NAME};
#my $cas_ticket;
#$cas_ticket = $FORM->param('ticket');
#$USER = undef;
#($USER) = CoGe::Accessory::Web->login_cas(
#    ticket   => $cas_ticket,
#    coge     => $coge,
#    this_url => $FORM->url()
#) if ($cas_ticket);
#($USER) = CoGe::Accessory::Web->get_user(
#    cookie_name => $COOKIE_NAME,
#    coge        => $coge
#) unless $USER;

my $dsgid = $FORM->param('dsgid');
my $ftid  = $FORM->param('ftid');
my $prot  = $FORM->param('p');
unless ( $dsgid ) {
    print $FORM->header;
    print "No genome id specified.\n";
    exit;
}
unless ($ftid) {
    print $FORM->header;
    print "No feature type id specified.\n";
    exit;
}

my $dsg = $coge->resultset('Genome')->find($dsgid);

if ( !$USER->has_access_to_genome($dsg) ) {
    print $FORM->header;
    print $USER->name,"\n";    
    print "Permission denied";
    exit;
}

my $ft = $coge->resultset('FeatureType')->find($ftid);
my $file_name;
$file_name .= $dsgid;
$file_name .= "-" . $ft->name;
$file_name .= "-prot" if $prot;
$file_name .= ".fasta";

print qq{Content-Type: application/force-download
Content-Disposition: attachment; filename="$file_name"

};

my $count = 1;
my @feats = $dsg->features( { feature_type_id => $ftid } );
my $org = $dsg->organism->name;

foreach my $feat (@feats) {
    my ($chr) = $feat->chromosome;    #=~/(\d+)/;
    my $name;

    eval {{
        foreach my $n ( $feat->names ) {
            $name = $n;
            last unless $name =~ /\s/;
        }
        $name =~ s/\s+/_/g;

        my $title = join( "||",
            $org,              $chr,      $feat->start,
            $feat->stop,       $name,     $feat->strand,
            $feat->type->name, $feat->id, $count );
        $title = ">" . $title;

        if ($prot) {
            my (@seqs) = $feat->protein_sequence( dsgid => $dsg->id );
            unless (scalar @seqs) {
                my $msg = "Skipping feature $name (id=" . $feat->id . ") due to no seqs";
                warn $msg;
                print '# ', $msg, "\n";
                next;
            }

            for (my $i = 0; $i < @seqs; $i++) {
                print $title, "||frame$i\n", $seqs[$i], "\n";
            }
        }
        else {
            print $title, "\n", $feat->genomic_sequence, "\n";
        }

        $count++;
    }};

    if ($@) {
        say STDERR "DOWNLOAD ABORTED: $file_name";
        last;
    }
}
