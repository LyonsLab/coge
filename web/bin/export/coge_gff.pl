#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use CGI;

no warnings 'redefine';

my $form = new CGI;

my ( $coge, $user ) = CoGe::Accessory::Web->init( cgi => $form );

my $dsgid   = $form->param('dsgid');
my $dsid    = $form->param('dsid');
my $id_type = $form->param('id_type');
my $annos   = 0;                         #flag for printing annotations
$annos = $form->param('annos') if $form->param('annos');
my $cds = 0;    #flag for printing only genes, mRNA, and CDSs
$cds = $form->param('cds') if $form->param('cds');
my $name_unique = 0;
$name_unique = $form->param('nu') if $form->param('nu');
my $upa = $form->param('upa') if $form->param('upa'); #unqiue_parent_annotations

my $item;
my $id;
my $org;
if ($dsgid) {
    $item = $coge->resultset('Genome')->find($dsgid);
    if ( $item && $item->restricted ) {
        if ( !$user->has_access_to_genome($item) ) {
            print $form->header;
            print "Error: permission denied\n";
            exit;
        }
    }
    $id  = $dsgid;
    $org = $item->organism->name . "_dsgid";
}
elsif ($dsid) {
    $item = $coge->resultset('Dataset')->find($dsid);
    $id   = $dsid;
    $org  = $item->organism->name . "_dsid";
}
unless ($item) {
    print $form->header('text');
    print "Unable to complete database transaction.\n";
    print
"Please contact coge administrator with a description of your problem for additional help.\n";
    exit();
}

$org =~ s/\///g;
$org =~ s/\s+/_/g;
$org =~ s/\(//g;
$org =~ s/\)//g;
$org =~ s/://g;
$org =~ s/;//g;
$org =~ s/#/_/g;
$org =~ s/'//g;
$org =~ s/"//g;

print "Content-disposition: attachement; filename=", $org, $id, ".gff\n\n";
$item->gff(
    print                     => 1,
    annos                     => $annos,
    cds                       => $cds,
    name_unique               => $name_unique,
    id_type                   => $id_type,
    unique_parent_annotations => $upa
);
