#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use CoGe::Accessory::Web;
use Getopt::Long;

use vars qw($dsgid $conffile);
GetOptions(
	   "dsgid=i"=>\$dsgid,
	   "config_file|cf=s"=>\$conffile,
	  );

unless ( $dsgid )
  {
    print qq{Usage:  $0 <dsgid>

Will generate a dump of all features for a dataset_group in bed format
};
    exit;
  }

my $P = CoGe::Accessory::Web::get_defaults($conffile);
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );
my $dsg = $coge->resultset('Genome')->find($dsgid);
foreach my $ds ($dsg->datasets)
  {
    my $rs = $coge->resultset('Feature');
    $rs->result_class('DBIx::Class::ResultClass::HashRefInflator');
    foreach my $ref ($rs->search({dataset_id => $ds->id, feature_type_id=>[3,5,8]}))
      {
	print join ("\t", $ref->chromosome, $ref->start, $ref->stop, $ref->feature_id),"\n";
      }
  }
