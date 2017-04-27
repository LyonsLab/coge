#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Sort::Versions;
use CoGe::Accessory::Web;
use CGI;

no warnings 'redefine';

my $P      = CoGe::Accessory::Web::get_defaults();
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};

my $connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
my $coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

my $COOKIE_NAME  = $P->{COOKIE_NAME};
my $FORM         = new CGI;
my ($cas_ticket) = $FORM->param('ticket');
my $USER         = undef;
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

my $dsgid = $FORM->param('dsgid');

my $dsg = $coge->resultset('Genome')->find($dsgid);
unless ($dsg) {
    print $FORM->header('text');
    print "Unable to complete database transaction.\n";
    print
"Please contact coge administrator with a description of your problem for additional help.\n";
    exit();
}

if ( !$USER->has_access_to_genome($dsg) ) {
    print $FORM->header('text');
    print
"Unauthorized access to restricted data.  Please contact coge administrator for additional help.\n";
    exit();
}
my $org = $dsg->organism->name;
$org =~ s/\s+/_/g;
my $header = "Content-disposition: attachement; filename=";    #test.gff\n\n";
$header .= $org;
$header .= "dsgid$dsgid";
$header .= "_tbl.txt\n\n";
print $header;

my %chr2ds;
foreach my $ds ( $dsg->datasets ) {
    map { $chr2ds{$_} = $ds } $ds->chromosomes;
}
foreach my $chr ( sort { versioncmp( $a, $b ) } $dsg->chromosomes ) {
    print ">Features $chr\n";
    foreach my $feat (
        sort {
            $a->start <=> $b->start
              || $a->feature_type_id <=> $b->feature_type_id
        } $chr2ds{$chr}
        ->features( { chromosome => $chr }, { "order_by" => "start ASC" } )
      )
    {
        next if $feat->type->name eq "chromosome";
        next if $feat->type->name =~ /utr/i;

        my $locs = get_locs($feat);

        my $item = shift @$locs;
        print join( "\t", @$item, $feat->type->name ), "\n";
        foreach my $loc (@$locs) {
            print join( "\t", @$loc ), "\n";
        }
        my $name_tag = get_name_tag($feat);
        my ( $pri_name, @names ) = $feat->names;
        print "\t\t\t", $name_tag, "\t";
        if ( $name_tag =~ /_id/ ) {
            print "gnl|dbname|";
        }
        print $pri_name, "\n";

        foreach my $name (@names) {
            print "\t\t\t", join( "\t", "alt_name", $name ), "\n";
        }
        foreach my $anno ( $feat->annotations ) {
            print "\t\t\t", join( "\t", $anno->type->name, $anno->annotation ),
              "\n";
        }
    }
}

sub get_name_tag {
    my $feat = shift;
    my $tag;
    if    ( $feat->type->name =~ /mRNA/ ) { $tag = "transcript_id" }
    elsif ( $feat->type->name =~ /CDS/ )  { $tag = "protein_id" }
    else                                  { $tag = "locus_tag" }
    return $tag;
}

sub get_locs {
    my $feat = shift;
    my @locs =
      map { [ $_->start, $_->stop ] }
      sort { $a->start <=> $b->start } $feat->locations;
    if ( $feat->strand =~ /-/ ) {
        @locs = reverse @locs;
        foreach my $item (@locs) {
            ( $item->[0], $item->[1] ) = ( $item->[1], $item->[0] );
        }
    }
    return \@locs;
}
