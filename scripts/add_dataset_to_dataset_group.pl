#!/usr/bin/perl -w
use strict;
use CoGeX;
use Getopt::Long;

my $dsid;
my $dsgid;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

GetOptions ( "dsid=i"=>\$dsid,
	     "dsgid=i"=>\$dsgid,
	     );
unless ($dsid && $dsgid)
  {
    print qq{
USAGE  $0 -dsid <dataset id> -dsgid <dataset group id>

Purpose:  Add the approprate row to the dataset_connector table so a dataset is in a dataset_group
};
    exit;
  }

my $res = $coge->resultset('DatasetConnector')->find_or_create({dataset_id=>$dsid, dataset_group_id=>$dsgid});

if ($res)
  {
    print "Created dataset connector with database id: ".$res->id,"\n";
  }
else
  {
    print "Failed to create dataset connector\n";
  }
