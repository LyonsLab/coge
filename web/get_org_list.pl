#!/usr/bin/perl -w

use strict;
use CoGeX;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::Web;
no warnings 'redefine';

my $P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};

my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr =
  "dbi:$P->{DB}:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
my $coge  = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
my $FORM  = new CGI;
my $dsgid = $FORM->param('dsgid');
my $oid   = $FORM->param('oid');
my $name  = $FORM->param('name');
my $desc  = $FORM->param('desc');

my @orgs;

if ( $name || $desc ) {
    my $search = {};
    $search->{description} = { like => "%" . $desc . "%" } if $desc;
    $search->{name}        = { like => "%" . $name . "%" } if $name;
    @orgs = $coge->resultset("Organism")->search($search);
}
elsif ($oid) {
    my @orgs = $coge->resultset("Organism")->find($oid);
}
elsif ($dsgid) {
    my $dsg = $coge->resultset("Genome")->find($dsgid);
    @orgs = $dsg->organism;
}
my $header =
  "Content-disposition: attachement; filename=CoGe_Organism_List.txt\n\n";
print $header;

#print $FORM->header('text');
#print "Context-Type:text/plain\n\n";
foreach my $org ( sort { $a->name cmp $b->name } @orgs ) {
    print join( "\t", $org->name, $org->description ), "\n";
}
