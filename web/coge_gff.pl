#!/usr/bin/perl -w

use CoGeX;
use CoGe::Accessory::Web;
use strict;
use CGI;
use CoGe::Accessory::LogUser;

no warnings 'redefine';
my $P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );
my $FORM = new CGI;
my $dsgid = $FORM->param('dsgid');
my $dsid = $FORM->param('dsid');
my $id_type = $FORM->param('id_type');
my $annos = 0; #flag for printing annotations
$annos = $FORM->param('annos') if $FORM->param('annos');
my $cds = 0; #flag for printing only genes, mRNA, and CDSs
$cds = $FORM->param('cds') if $FORM->param('cds');
my $name_unique = 0;
$name_unique = $FORM->param('nu') if $FORM->param('nu');

my ($USER) = CoGe::Accessory::LogUser->get_user();

my $item;
my $id;
my $org;
if ($dsgid)
  {
    $item = $coge->resultset('DatasetGroup')->find($dsgid);
    next if $USER->user_name =~ /public/i && $item->organism->restricted;
    $id = $dsgid;
    $org = $item->organism->name."_dsgid";
  }
elsif ($dsid)
  {
    $item = $coge->resultset('Dataset')->find($dsid);
    next if $USER->user_name =~ /public/i && $item->organism->restricted;
    $id = $dsid;
    $org = $item->organism->name."_dsid";
  }
unless ($item)
  {
    print $FORM->header('text');
    print "Unable to complete database transaction.\n";
    print "Please contact coge administrator with a description of your problem for additional help.\n";
    exit();
  }


$org =~ s/\s+/_/g;
my $header = "Content-disposition: attachement; filename=";#test.gff\n\n";
$header .= $org;
$header .= $id;
$header .=".gff\n\n";
print $header;
$item->gff(print=>1, annos=>$annos, cds=>$cds, name_unique=>$name_unique, id_type=>$id_type);
